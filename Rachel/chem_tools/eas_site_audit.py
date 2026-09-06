"""Conservative, knowledge-backed EAS site-selectivity audit."""

from __future__ import annotations

import re
from typing import Any, Dict, Mapping, Optional, Sequence

from rdkit import Chem

from ._rdkit_utils import parse_mol


_SUPPORTED_EAS_TERMS = ("aromatic nitration", "arene nitration", "nitration")
_ARYL_NITRO_QUERY = Chem.MolFromSmarts("[c:1]-[N+:2](=[O:3])[O-:4]")

_POSITION_BY_DISTANCE = {1: "ortho", 2: "meta", 3: "para"}
_REQUIRED_EVIDENCE = [
    "substrate-specific literature precedent, experimental result, or internal SOP",
    "conditions and observed regioisomer distribution for a close structural analogue",
]


def _is_eas_action(reaction_type: str, action_context: Optional[Mapping[str, Any]]) -> bool:
    fields = [reaction_type]
    if isinstance(action_context, Mapping):
        fields.extend(
            str(action_context.get(key, "") or "")
            for key in ("reaction_type", "reaction_name", "template_name")
        )
    normalized = " ".join(fields).lower().replace("_", " ")
    normalized = re.sub(r"\s+", " ", normalized)
    return any(
        re.search(rf"(?<!\w){re.escape(term)}(?!\w)", normalized)
        for term in _SUPPORTED_EAS_TERMS
    )


def _not_applicable(reason_code: str, summary: str) -> Dict[str, Any]:
    return {
        "ok": True,
        "applicable": False,
        "state": "not_applicable",
        "reason_code": reason_code,
        "summary": summary,
    }


def _inconclusive(reason_code: str, summary: str) -> Dict[str, Any]:
    return {
        "ok": True,
        "applicable": True,
        "state": "inconclusive",
        "reason_code": reason_code,
        "method": "heuristic_director_consensus",
        "atom_mapping_is_selectivity_proof": False,
        "required_evidence": list(_REQUIRED_EVIDENCE),
        "summary": summary,
    }


def _ring_distance(ring: Sequence[int], first: int, second: int) -> int:
    first_pos = ring.index(first)
    second_pos = ring.index(second)
    direct = abs(first_pos - second_pos)
    return min(direct, len(ring) - direct)


def _query_map_index(query: Chem.Mol, atom_map: int) -> Optional[int]:
    for atom in query.GetAtoms():
        if atom.GetAtomMapNum() == atom_map:
            return atom.GetIdx()
    return None


def audit_eas_site_selectivity(
    product_smiles: str,
    precursors: Sequence[str],
    *,
    reaction_type: str = "",
    action_context: Optional[Mapping[str, Any]] = None,
    site_audit: Optional[Mapping[str, Any]] = None,
    knowledge_profile=None,
) -> Dict[str, Any]:
    """Assess director support without treating the heuristic as selectivity proof."""
    if not _is_eas_action(reaction_type, action_context):
        return _not_applicable(
            "reaction_family_not_eas",
            "The initial EAS site audit applies only to explicitly declared nitration.",
        )

    site_audit = site_audit if isinstance(site_audit, Mapping) else {}
    if not site_audit.get("site_retentive"):
        return _inconclusive(
            "same_core_site_map_unavailable",
            "EAS was declared, but no same-core aromatic site map was available.",
        )
    changed_sites = site_audit.get("changed_sites") or []
    if len(changed_sites) != 1 or int(site_audit.get("changed_site_count", -1)) != 1:
        return _inconclusive(
            "single_substitution_site_unavailable",
            "EAS site selectivity requires one unambiguous changed aromatic site.",
        )
    product_mol = parse_mol(product_smiles)
    if product_mol is None or _ARYL_NITRO_QUERY is None:
        return _inconclusive(
            "desired_site_atom_unavailable",
            "The nitration product could not be assigned to a molecular graph.",
        )
    product_nitro_matches = product_mol.GetSubstructMatches(
        _ARYL_NITRO_QUERY, uniquify=True
    )
    precursor_nitro_count = 0
    for precursor in precursors:
        precursor_mol = parse_mol(precursor)
        if precursor_mol is not None:
            precursor_nitro_count += len(
                precursor_mol.GetSubstructMatches(_ARYL_NITRO_QUERY, uniquify=True)
            )
    if len(product_nitro_matches) != 1 or precursor_nitro_count != 0:
        return _inconclusive(
            "unique_new_aryl_nitro_group_unavailable",
            "The initial nitration audit requires exactly one product aryl nitro group "
            "and no aryl nitro group in the precursor set.",
        )
    desired_atom_idx = int(product_nitro_matches[0][0])

    candidate_rings = []
    for atom_ring in product_mol.GetRingInfo().AtomRings():
        ring = list(atom_ring)
        if desired_atom_idx not in ring or len(ring) != 6:
            continue
        if all(
            product_mol.GetAtomWithIdx(atom_idx).GetIsAromatic()
            and product_mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6
            for atom_idx in ring
        ):
            candidate_rings.append(ring)
    if len(candidate_rings) != 1:
        return _inconclusive(
            "supported_aromatic_ring_unavailable",
            "The initial EAS audit supports only one unambiguous six-membered carbocycle.",
        )
    ring = candidate_rings[0]

    ring_atom_ids = set(ring)
    preserved_atom_ids = {
        atom_idx
        for atom_idx in ring
        if atom_idx != desired_atom_idx
        and any(
            neighbor.GetIdx() not in ring_atom_ids
            for neighbor in product_mol.GetAtomWithIdx(atom_idx).GetNeighbors()
        )
    }
    available_sites = {
        atom_idx
        for atom_idx in ring
        if atom_idx == desired_atom_idx
        or product_mol.GetAtomWithIdx(atom_idx).GetTotalNumHs() > 0
    }
    if not preserved_atom_ids:
        return _not_applicable(
            "no_preserved_aromatic_substituents",
            "The audited ring has no preserved substituent, so a director-based "
            "regioselectivity audit is not applicable.",
        )

    if knowledge_profile is None:
        from Rachel.knowledge import get_base_profile

        knowledge_profile = get_base_profile()
    resource = knowledge_profile.get("chem.eas_directing_rules")
    entries = resource.get("entries") if isinstance(resource, Mapping) else None
    if not isinstance(entries, list):
        return _inconclusive(
            "eas_rule_resource_unavailable",
            "The EAS directing-rule resource is unavailable.",
        )

    directors = []
    knowledge_refs = []
    seen_directors = set()
    for entry in sorted(
        (item for item in entries if isinstance(item, Mapping)),
        key=lambda item: (-int(item.get("priority", 0) or 0), str(item.get("id", ""))),
    ):
        entry_id = str(entry.get("id", "") or "")
        query = Chem.MolFromSmarts(str(entry.get("smarts", "") or ""))
        if query is None:
            return _inconclusive(
                "invalid_eas_directing_rule",
                f"EAS directing rule could not be parsed: {entry_id}",
            )
        try:
            ring_atom_map = int(entry.get("ring_atom_map", 0))
        except (TypeError, ValueError):
            ring_atom_map = 0
        query_ring_idx = _query_map_index(query, ring_atom_map)
        if query_ring_idx is None:
            return _inconclusive(
                "invalid_eas_directing_rule",
                f"EAS directing rule lacks its mapped ring atom: {entry_id}",
            )
        for match in product_mol.GetSubstructMatches(query, uniquify=True):
            ring_atom_idx = int(match[query_ring_idx])
            key = (entry_id, ring_atom_idx)
            if (
                ring_atom_idx not in ring
                or ring_atom_idx not in preserved_atom_ids
                or key in seen_directors
            ):
                continue
            seen_directors.add(key)
            knowledge_ref = knowledge_profile.source(
                "chem.eas_directing_rules", entry_id
            )
            if knowledge_ref not in knowledge_refs:
                knowledge_refs.append(knowledge_ref)
            directors.append(
                {
                    "rule_id": entry_id,
                    "name": str(entry.get("name", "") or entry_id),
                    "ring_atom_idx": ring_atom_idx,
                    "directing_positions": list(entry.get("directing_positions", []) or []),
                    "electronic_effect": str(entry.get("electronic_effect", "") or ""),
                    "knowledge_ref": knowledge_ref,
                }
            )

    if not directors:
        result = _inconclusive(
            "no_recognized_preserved_director",
            "No preserved substituent was covered by the active EAS directing rules.",
        )
        result["desired_atom_idx"] = desired_atom_idx
        return result

    covered_director_sites = {director["ring_atom_idx"] for director in directors}
    uncovered_preserved_sites = sorted(
        (preserved_atom_ids & set(ring)) - covered_director_sites
    )
    if uncovered_preserved_sites:
        result = _inconclusive(
            "preserved_substituent_rule_coverage_incomplete",
            "At least one preserved aromatic substituent is outside the active EAS "
            "directing-rule coverage, so no site verdict was issued.",
        )
        result.update(
            {
                "desired_atom_idx": desired_atom_idx,
                "uncovered_preserved_atom_ids": uncovered_preserved_sites,
                "directors": directors,
                "knowledge_refs": knowledge_refs,
            }
        )
        return result

    site_rows = []
    for site_atom_idx in sorted(available_sites):
        supporting_rule_ids = []
        for director in directors:
            distance = _ring_distance(ring, director["ring_atom_idx"], site_atom_idx)
            position = _POSITION_BY_DISTANCE.get(distance, "")
            if position in director["directing_positions"]:
                supporting_rule_ids.append(director["rule_id"])
        site_rows.append(
            {
                "atom_idx": site_atom_idx,
                "is_desired": site_atom_idx == desired_atom_idx,
                "supporting_rule_ids": supporting_rule_ids,
                "support_count": len(supporting_rule_ids),
            }
        )

    desired_row = next(row for row in site_rows if row["is_desired"])
    alternative_rows = [row for row in site_rows if not row["is_desired"]]
    supported_alternatives = [
        row for row in alternative_rows if row["support_count"] > 0
    ]
    if desired_row["support_count"] == 0:
        state = "unsupported"
        reason_code = "desired_site_not_supported_by_director_consensus"
        summary = (
            "The desired EAS site is not supported by the active director consensus; "
            "substrate-specific evidence is required."
        )
    elif supported_alternatives:
        state = "competing"
        reason_code = "supported_site_has_competition"
        summary = (
            "The desired EAS site has director support, but other aromatic C-H sites are "
            "also supported and may compete."
        )
    else:
        state = "supported"
        reason_code = "desired_site_supported_without_detected_competition"
        summary = (
            "The desired EAS site is supported by the active director rules; this remains "
            "heuristic evidence, not experimental selectivity proof."
        )

    result = {
        "ok": True,
        "applicable": True,
        "state": state,
        "reason_code": reason_code,
        "method": "heuristic_director_consensus",
        "atom_mapping_is_selectivity_proof": False,
        "desired_atom_idx": desired_atom_idx,
        "directors": directors,
        "candidate_sites": site_rows,
        "knowledge_refs": knowledge_refs,
        "summary": summary,
    }
    if state == "unsupported":
        result["required_evidence"] = list(_REQUIRED_EVIDENCE)
    return result
