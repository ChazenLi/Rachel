"""Validation and lifecycle management for non-executable knowledge packs."""

from __future__ import annotations

import hashlib
import json
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping

from .conditions import KnowledgeConditionError, condition_matches
from .profile import RESOURCE_SCHEMA, SCHEMA


_PACK_PART_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_BASE_MANIFEST = Path(__file__).resolve().parent / "base" / "manifest.json"
_ROUTE_EVIDENCE_KINDS = {
    "route_review",
    "decision_audit",
    "experimental_outcome",
}
_HIGH_EVIDENCE_RESOURCES = {
    "chem.reactions",
    "chem.functional_groups",
    "chem.reactive_sites",
    "chem.protecting_groups",
    "chem.smart_cap_rules",
    "chem.reaction_families",
    "chem.eas_directing_rules",
    "chem.cs_smarts",
    "chem.dangerous_combos",
    "chem.fg_compatibility",
    "chem.selectivity_conflicts",
    "chem.selectivity_reactivity",
    "chem.structural_alerts",
}
_AUTHORITATIVE_EVIDENCE_KINDS = {
    "experimental_run",
    "internal_run",
    "literature",
    "internal_sop",
    "sop",
}


class KnowledgePackError(ValueError):
    """Raised when a pack-management command cannot complete safely."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_json(path: Path, *, label: str) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise KnowledgePackError(f"invalid {label}: {path.name}") from exc


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = ""
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_name = stream.name
            json.dump(payload, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        Path(temporary_name).replace(path)
    finally:
        if temporary_name:
            temporary = Path(temporary_name)
            if temporary.exists():
                temporary.unlink()


def _manifest_path(pack: Path | str) -> Path:
    path = Path(pack).resolve()
    if path.is_dir():
        path = path / "manifest.json"
    if not path.is_file():
        raise KnowledgePackError(f"knowledge pack manifest not found: {path}")
    return path


def _validate_smarts_fields(value: Any, *, resource: str, entry_id: str) -> None:
    try:
        from rdkit import Chem
    except ImportError as exc:  # pragma: no cover - Rachel runtime requires RDKit
        raise KnowledgePackError("RDKit is required to validate chemistry resources") from exc

    if isinstance(value, dict):
        for key, child in value.items():
            if key == "smarts" or key.endswith("_smarts") or key == "prereq":
                if not isinstance(child, str) or Chem.MolFromSmarts(child) is None:
                    raise KnowledgePackError(
                        f"invalid SMARTS: {resource}/{entry_id}/{key}"
                    )
            else:
                _validate_smarts_fields(child, resource=resource, entry_id=entry_id)
    elif isinstance(value, list):
        for child in value:
            _validate_smarts_fields(child, resource=resource, entry_id=entry_id)


def _validate_reaction_entry(entry: Mapping[str, Any], *, entry_id: str) -> None:
    value = entry.get("value")
    required = ("name", "type", "prereq", "rxn")
    if not isinstance(value, dict) or any(
        not isinstance(value.get(field), str) or not value.get(field).strip()
        for field in required
    ):
        raise KnowledgePackError(
            f"reaction entry requires name, type, prereq, and rxn: {entry_id}"
        )
    if value["type"] not in {"forward", "retro"}:
        raise KnowledgePackError(f"reaction entry type must be forward or retro: {entry_id}")
    try:
        from rdkit.Chem import rdChemReactions
    except ImportError as exc:  # pragma: no cover - Rachel runtime requires RDKit
        raise KnowledgePackError("RDKit is required to validate reaction resources") from exc
    try:
        reaction = rdChemReactions.ReactionFromSmarts(value["rxn"])
    except (RuntimeError, ValueError):
        reaction = None
    if reaction is None:
        raise KnowledgePackError(f"invalid reaction SMARTS: {entry_id}")


def _validate_eas_directing_rule(
    entry: Mapping[str, Any], *, entry_id: str
) -> None:
    from rdkit import Chem

    name = entry.get("name")
    smarts = entry.get("smarts")
    if not isinstance(name, str) or not name.strip() or not isinstance(smarts, str):
        raise KnowledgePackError(
            f"EAS directing rule requires name and smarts: {entry_id}"
        )
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        raise KnowledgePackError(f"invalid EAS directing SMARTS: {entry_id}")
    ring_atom_map = entry.get("ring_atom_map")
    if isinstance(ring_atom_map, bool) or not isinstance(ring_atom_map, int):
        raise KnowledgePackError(
            f"EAS directing rule requires integer ring_atom_map: {entry_id}"
        )
    mapped_atoms = [
        atom for atom in query.GetAtoms() if atom.GetAtomMapNum() == ring_atom_map
    ]
    if (
        len(mapped_atoms) != 1
        or not mapped_atoms[0].GetIsAromatic()
        or mapped_atoms[0].GetAtomicNum() != 6
    ):
        raise KnowledgePackError(
            f"EAS ring_atom_map must identify one aromatic carbon: {entry_id}"
        )
    positions = entry.get("directing_positions")
    if (
        not isinstance(positions, list)
        or not positions
        or len(set(positions)) != len(positions)
        or not all(position in {"ortho", "meta", "para"} for position in positions)
    ):
        raise KnowledgePackError(f"invalid EAS directing positions: {entry_id}")
    if entry.get("electronic_effect") not in {
        "activating",
        "deactivating",
        "context_dependent",
    }:
        raise KnowledgePackError(f"invalid EAS electronic effect: {entry_id}")
    priority = entry.get("priority")
    if isinstance(priority, bool) or not isinstance(priority, int):
        raise KnowledgePackError(f"invalid EAS rule priority: {entry_id}")


def _validate_conditions(entry: Mapping[str, Any], *, entry_id: str) -> None:
    for field in ("condition", "activation"):
        if field not in entry:
            continue
        try:
            condition_matches(entry[field], {})
        except (KnowledgeConditionError, TypeError, ValueError) as exc:
            raise KnowledgePackError(
                f"unsupported declarative condition: {entry_id}/{field}"
            ) from exc


def _declares_relationship(entry: Mapping[str, Any], other_id: str) -> bool:
    metadata = entry.get("_knowledge")
    if not isinstance(metadata, dict):
        return False
    for field in ("coexist", "supersedes"):
        value = metadata.get(field)
        if isinstance(value, str) and value == other_id:
            return True
        if isinstance(value, list) and other_id in value:
            return True
    return False


def _relationship_resolves(
    first: Mapping[str, Any], second: Mapping[str, Any]
) -> bool:
    first_id = str(first.get("id", ""))
    second_id = str(second.get("id", ""))
    return _declares_relationship(first, second_id) or _declares_relationship(
        second, first_id
    )


def _prompt_effect_conflicts(entries: Mapping[str, Any]) -> list[Dict[str, Any]]:
    positive = {"allow", "enable", "require", "encourage", "prefer"}
    negative = {"deny", "disable", "forbid", "discourage", "avoid"}
    scoped: Dict[str, list[Mapping[str, Any]]] = {}
    for entry in entries.values():
        if not isinstance(entry, dict) or "scope" not in entry:
            continue
        effect = str(entry.get("effect", "")).lower()
        if effect not in positive | negative:
            continue
        scope = json.dumps(
            entry["scope"], sort_keys=True, ensure_ascii=True, separators=(",", ":")
        )
        scoped.setdefault(scope, []).append(entry)

    conflicts: list[Dict[str, Any]] = []
    for group in scoped.values():
        for index, first in enumerate(group):
            first_effect = str(first["effect"]).lower()
            for second in group[index + 1 :]:
                second_effect = str(second["effect"]).lower()
                opposite = (first_effect in positive and second_effect in negative) or (
                    first_effect in negative and second_effect in positive
                )
                if not opposite or _relationship_resolves(first, second):
                    continue
                conflicts.append(
                    {
                        "type": "opposite_prompt_effect",
                        "resource": "prompt.runtime_policy",
                        "entry_ids": sorted([str(first["id"]), str(second["id"])]),
                    }
                )
    return conflicts


def _risk_conflicts(
    resources: Mapping[str, Mapping[str, Any]]
) -> list[Dict[str, Any]]:
    risk_resources = {
        "chem.dangerous_combos",
        "chem.fg_compatibility",
        "chem.selectivity_conflicts",
        "chem.selectivity_reactivity",
        "chem.structural_alerts",
    }
    conflicts: list[Dict[str, Any]] = []
    for resource in sorted(risk_resources):
        by_key: Dict[str, list[Mapping[str, Any]]] = {}
        for entry in resources.get(resource, {}).values():
            if not isinstance(entry, dict):
                continue
            risk_key = str(entry.get("risk_key", "")).strip()
            severity = str(entry.get("severity", "")).strip().lower()
            if risk_key and severity:
                by_key.setdefault(risk_key, []).append(entry)
        for group in by_key.values():
            for index, first in enumerate(group):
                for second in group[index + 1 :]:
                    if first["severity"] == second["severity"]:
                        continue
                    if _relationship_resolves(first, second):
                        continue
                    conflicts.append(
                        {
                            "type": "risk_severity_mismatch",
                            "resource": resource,
                            "entry_ids": sorted(
                                [str(first["id"]), str(second["id"])]
                            ),
                        }
                    )
    return conflicts


def _family_conflicts(entries: Mapping[str, Any]) -> list[Dict[str, Any]]:
    by_family: Dict[str, list[Mapping[str, Any]]] = {}
    for entry in entries.values():
        if not isinstance(entry, dict):
            continue
        family_key = str(entry.get("family_key", "")).strip()
        if family_key:
            by_family.setdefault(family_key, []).append(entry)

    conflicts: list[Dict[str, Any]] = []
    for group in by_family.values():
        for entry in group:
            allowed = set(entry.get("allowed_deltas", []) or [])
            forbidden = set(entry.get("forbidden_deltas", []) or [])
            if allowed & forbidden:
                conflicts.append(
                    {
                        "type": "family_delta_disagreement",
                        "resource": "chem.reaction_families",
                        "entry_ids": [str(entry["id"])],
                    }
                )
        for index, first in enumerate(group):
            first_allowed = set(first.get("allowed_deltas", []) or [])
            first_forbidden = set(first.get("forbidden_deltas", []) or [])
            for second in group[index + 1 :]:
                second_allowed = set(second.get("allowed_deltas", []) or [])
                second_forbidden = set(second.get("forbidden_deltas", []) or [])
                disagreement = (first_allowed & second_forbidden) or (
                    second_allowed & first_forbidden
                )
                if not disagreement or _relationship_resolves(first, second):
                    continue
                conflicts.append(
                    {
                        "type": "family_delta_disagreement",
                        "resource": "chem.reaction_families",
                        "entry_ids": sorted([str(first["id"]), str(second["id"])]),
                    }
                )
    return conflicts


class KnowledgePackManager:
    """Deep module for validating and publishing immutable knowledge packs."""

    def validate(self, pack: Path | str) -> Dict[str, Any]:
        """Validate one published external pack without activating it."""
        manifest_path = _manifest_path(pack)
        manifest = _read_json(manifest_path, label="knowledge manifest")
        if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
            raise KnowledgePackError("unsupported knowledge manifest schema")

        pack_id = str(manifest.get("id", ""))
        version = str(manifest.get("version", ""))
        kind = str(manifest.get("kind", ""))
        if not _PACK_PART_RE.fullmatch(pack_id):
            raise KnowledgePackError("knowledge pack id is invalid")
        if not _PACK_PART_RE.fullmatch(version):
            raise KnowledgePackError("knowledge pack version is invalid")
        if kind not in {"team", "project"}:
            raise KnowledgePackError("external knowledge pack kind must be team or project")

        dependencies = manifest.get("dependencies")
        if not isinstance(dependencies, list) or not all(
            isinstance(item, str) for item in dependencies
        ):
            raise KnowledgePackError("knowledge pack dependencies must be versioned references")
        for dependency in dependencies:
            self._split_pack_ref(dependency)

        resources = manifest.get("resources")
        if not isinstance(resources, dict) or not resources:
            raise KnowledgePackError("knowledge pack resources must be a non-empty object")
        base_manifest = _read_json(_BASE_MANIFEST, label="base knowledge manifest")
        base_resources = base_manifest.get("resources", {})

        entry_count = 0
        pack_dir = manifest_path.parent.resolve()
        for resource, spec in resources.items():
            if resource not in base_resources:
                raise KnowledgePackError(f"unknown knowledge resource: {resource}")
            if not isinstance(spec, Mapping):
                raise KnowledgePackError(f"invalid resource declaration: {resource}")
            resource_path = (pack_dir / str(spec.get("path", ""))).resolve()
            try:
                resource_path.relative_to(pack_dir)
            except ValueError as exc:
                raise KnowledgePackError(
                    f"knowledge resource escapes its pack directory: {resource}"
                ) from exc
            if resource_path.suffix.lower() != ".json" or not resource_path.is_file():
                raise KnowledgePackError(
                    f"knowledge resource must be a JSON file: {resource}"
                )
            expected_hash = str(spec.get("sha256", "")).lower()
            if not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
                raise KnowledgePackError(f"invalid resource sha256: {resource}")
            if _sha256(resource_path) != expected_hash:
                raise KnowledgePackError(f"knowledge resource hash mismatch: {resource}")

            overlay = _read_json(resource_path, label=f"knowledge resource {resource}")
            if not isinstance(overlay, dict) or overlay.get("schema") != RESOURCE_SCHEMA:
                raise KnowledgePackError(f"invalid knowledge overlay schema: {resource}")
            entries = overlay.get("entries")
            if not isinstance(entries, list):
                raise KnowledgePackError(f"knowledge overlay requires entries: {resource}")
            seen: set[str] = set()
            for entry in entries:
                if not isinstance(entry, dict):
                    raise KnowledgePackError(
                        f"knowledge overlay entry must be an object: {resource}"
                    )
                entry_id = str(entry.get("id", "")).strip()
                metadata = entry.get("_knowledge")
                operation = str(metadata.get("operation", "")) if isinstance(metadata, dict) else ""
                if not entry_id or operation not in {"add", "replace", "disable"}:
                    raise KnowledgePackError(
                        f"knowledge entry requires id and explicit operation: {resource}"
                    )
                if entry_id in seen:
                    raise KnowledgePackError(f"duplicate knowledge entry: {resource}/{entry_id}")
                seen.add(entry_id)
                if operation == "add" and not entry_id.startswith(
                    (f"{pack_id}.", f"{pack_id}:")
                ):
                    raise KnowledgePackError(
                        f"new knowledge entry must use pack namespace {pack_id}: {entry_id}"
                    )
                if operation in {"replace", "disable"}:
                    target = str(metadata.get("target", "")).strip()
                    supersedes = str(metadata.get("supersedes", "")).strip()
                    reason = str(metadata.get("reason", "")).strip()
                    if entry_id != target or supersedes != target or not reason:
                        raise KnowledgePackError(
                            f"{operation} requires target, matching supersedes, and reason: {resource}"
                        )
                if operation != "disable" and resource.startswith("chem."):
                    _validate_smarts_fields(
                        entry,
                        resource=resource,
                        entry_id=entry_id,
                    )
                if operation != "disable" and resource in {
                    "chem.smart_cap_rules",
                    "prompt.runtime_policy",
                }:
                    _validate_conditions(entry, entry_id=entry_id)
                if operation != "disable" and resource == "chem.reactions":
                    _validate_reaction_entry(entry, entry_id=entry_id)
                if operation != "disable" and resource == "chem.eas_directing_rules":
                    _validate_eas_directing_rule(entry, entry_id=entry_id)
                entry_count += 1

        return {
            "ok": True,
            "command": "validate",
            "pack": {"id": pack_id, "version": version, "kind": kind},
            "resource_count": len(resources),
            "entry_count": entry_count,
        }

    def diff(self, left: Path | str, right: Path | str) -> Dict[str, Any]:
        """Return a path-free, content-free entry diff for two valid packs."""
        left_summary = self.validate(left)
        right_summary = self.validate(right)
        left_entries = self._entries_by_resource(left)
        right_entries = self._entries_by_resource(right)

        resources: Dict[str, Dict[str, list[str]]] = {}
        totals = {"added": 0, "removed": 0, "changed": 0}
        for resource in sorted(set(left_entries) | set(right_entries)):
            before = left_entries.get(resource, {})
            after = right_entries.get(resource, {})
            added = sorted(set(after) - set(before))
            removed = sorted(set(before) - set(after))
            changed = sorted(
                entry_id
                for entry_id in set(before) & set(after)
                if before[entry_id] != after[entry_id]
            )
            if not (added or removed or changed):
                continue
            resources[resource] = {
                "added": added,
                "removed": removed,
                "changed": changed,
            }
            totals["added"] += len(added)
            totals["removed"] += len(removed)
            totals["changed"] += len(changed)

        return {
            "ok": True,
            "command": "diff",
            "left": left_summary["pack"],
            "right": right_summary["pack"],
            "summary": totals,
            "resources": resources,
        }

    def conflicts(self, packs: list[Path | str]) -> Dict[str, Any]:
        """Report blockers in an ordered base plus external-pack composition."""
        if not packs:
            raise KnowledgePackError("conflicts requires at least one pack")
        for pack in packs:
            self.validate(pack)

        base_manifest = _read_json(_BASE_MANIFEST, label="base knowledge manifest")
        base_specs = base_manifest["resources"]
        resolved = self._base_entries(base_manifest)
        conflicts: list[Dict[str, Any]] = []

        for pack in packs:
            for resource, entries in self._published_entries(pack).items():
                resource_entries = resolved.setdefault(resource, {})
                resource_locked = bool(base_specs[resource].get("locked"))
                for entry in entries:
                    entry_id = str(entry["id"])
                    metadata = entry["_knowledge"]
                    operation = str(metadata["operation"])
                    target = str(metadata.get("target", entry_id))
                    current = resource_entries.get(target)

                    if operation == "add":
                        if entry_id in resource_entries:
                            conflicts.append(
                                {
                                    "type": "duplicate_or_implicit_override",
                                    "resource": resource,
                                    "entry_ids": [entry_id],
                                }
                            )
                            continue
                        resource_entries[entry_id] = entry
                        continue

                    if current is None:
                        conflicts.append(
                            {
                                "type": "missing_operation_target",
                                "resource": resource,
                                "entry_ids": [target],
                            }
                        )
                        continue
                    if resource_locked:
                        conflicts.append(
                            {
                                "type": "locked_resource_modification",
                                "resource": resource,
                                "entry_ids": [target],
                            }
                        )
                        continue
                    if isinstance(current, dict) and current.get("locked") is True:
                        conflicts.append(
                            {
                                "type": "locked_entry_modification",
                                "resource": resource,
                                "entry_ids": [target],
                            }
                        )
                        continue
                    if operation == "disable":
                        resource_entries.pop(target)
                    else:
                        resource_entries[target] = entry

        conflicts.extend(
            _prompt_effect_conflicts(resolved.get("prompt.runtime_policy", {}))
        )
        conflicts.extend(_risk_conflicts(resolved))
        conflicts.extend(
            _family_conflicts(resolved.get("chem.reaction_families", {}))
        )

        conflicts.sort(
            key=lambda item: (
                item["resource"],
                item["type"],
                item["entry_ids"],
            )
        )
        return {
            "ok": not conflicts,
            "command": "conflicts",
            "pack_count": len(packs),
            "conflict_count": len(conflicts),
            "conflicts": conflicts,
        }

    def import_drafts(
        self,
        *,
        session: Path | str,
        workspace: Path | str,
        pack_id: str,
    ) -> Dict[str, Any]:
        """Import inactive session drafts into a path-free staging workspace."""
        if not _PACK_PART_RE.fullmatch(str(pack_id)):
            raise KnowledgePackError("target pack id is invalid")
        session_path = Path(session).resolve()
        session_payload = _read_json(session_path, label="Rachel session")
        if not isinstance(session_payload, dict):
            raise KnowledgePackError("Rachel session must contain a JSON object")
        session_id = str(session_payload.get("session_id", "")).strip()
        if not session_id or not _PACK_PART_RE.fullmatch(session_id):
            raise KnowledgePackError("Rachel session id is missing or invalid")
        audit_state = session_payload.get("audit_state")
        drafts = audit_state.get("knowledge_drafts") if isinstance(audit_state, dict) else None
        if not isinstance(drafts, list):
            raise KnowledgePackError("Rachel session has no knowledge drafts")
        tree = session_payload.get("tree")
        reaction_nodes = tree.get("reaction_nodes") if isinstance(tree, dict) else None
        if not isinstance(reaction_nodes, list):
            raise KnowledgePackError(
                "Rachel session route cannot verify knowledge draft sources"
            )
        route_actions_by_step: Dict[str, set[str]] = {}
        for reaction in reaction_nodes:
            if not isinstance(reaction, dict):
                continue
            step_id = str(reaction.get("step_id", "")).strip()
            if not step_id:
                continue
            action_ids: set[str] = set()
            llm_decision = reaction.get("llm_decision")
            if isinstance(llm_decision, dict):
                selected_candidate_id = str(
                    llm_decision.get("selected_candidate_id", "")
                ).strip()
                if selected_candidate_id:
                    action_ids.add(selected_candidate_id)
                decision_audit = llm_decision.get("decision_audit")
                if isinstance(decision_audit, dict):
                    selected_action_id = str(
                        decision_audit.get("selected_action_id", "")
                    ).strip()
                    if selected_action_id:
                        action_ids.add(selected_action_id)
                    rejected_actions = decision_audit.get("rejected_actions", [])
                    if isinstance(rejected_actions, list):
                        for rejected in rejected_actions:
                            if not isinstance(rejected, dict):
                                continue
                            rejected_action_id = str(
                                rejected.get("action_id", "")
                                or rejected.get("candidate_id", "")
                            ).strip()
                            if rejected_action_id:
                                action_ids.add(rejected_action_id)
            route_actions_by_step[step_id] = action_ids
        route_step_ids = set(route_actions_by_step)

        workspace_path = Path(workspace).resolve()
        if workspace_path.exists():
            staged = _read_json(workspace_path, label="knowledge staging workspace")
            if (
                not isinstance(staged, dict)
                or staged.get("schema") != "rachel.knowledge-staging.v1"
                or staged.get("pack_id") != pack_id
            ):
                raise KnowledgePackError("staging workspace identity mismatch")
        else:
            staged = {
                "schema": "rachel.knowledge-staging.v1",
                "pack_id": pack_id,
                "created_at": datetime.now(timezone.utc).isoformat(),
                "source_sessions": [],
                "drafts": [],
            }

        known_ids = {
            str(item.get("workspace_draft_id", ""))
            for item in staged.get("drafts", [])
            if isinstance(item, dict)
        }
        imported_ids: list[str] = []
        base_resources = _read_json(
            _BASE_MANIFEST, label="base knowledge manifest"
        )["resources"]
        for source in drafts:
            if not isinstance(source, dict) or source.get("target_pack_id") != pack_id:
                continue
            if source.get("status") != "staging" or source.get("active") is not False:
                raise KnowledgePackError("only inactive staging drafts can be imported")
            draft_id = str(source.get("draft_id", "")).strip()
            workspace_draft_id = f"{session_id}:{draft_id}"
            if not draft_id or workspace_draft_id in known_ids:
                continue
            resource = str(source.get("resource", ""))
            if resource not in base_resources:
                raise KnowledgePackError(f"unknown knowledge resource: {resource}")
            entry = source.get("entry")
            entry_id = str(entry.get("id", "")) if isinstance(entry, dict) else ""
            if not entry_id.startswith((f"{pack_id}.", f"{pack_id}:")):
                raise KnowledgePackError("draft entry id must use the target pack namespace")
            source_refs = source.get("source_refs")
            if not isinstance(source_refs, list) or not source_refs:
                raise KnowledgePackError("draft source_refs must be a non-empty list")
            for source_ref in source_refs:
                step_id = (
                    str(source_ref.get("step_id", "")).strip()
                    if isinstance(source_ref, dict)
                    else ""
                )
                if not step_id:
                    raise KnowledgePackError("draft source step is missing")
                if step_id not in route_step_ids:
                    raise KnowledgePackError(
                        f"draft source step not found: {step_id}"
                    )
                action_id = str(source_ref.get("action_id", "")).strip()
                if (
                    action_id
                    and action_id not in route_actions_by_step[step_id]
                ):
                    raise KnowledgePackError(
                        f"draft source action not found on step: {step_id}"
                    )
            imported = json.loads(json.dumps(source, ensure_ascii=False))
            imported["workspace_draft_id"] = workspace_draft_id
            imported["status"] = "staging"
            imported["active"] = False
            staged["drafts"].append(imported)
            known_ids.add(workspace_draft_id)
            imported_ids.append(workspace_draft_id)

        if not imported_ids:
            raise KnowledgePackError("no new staging drafts matched the target pack")
        knowledge = session_payload.get("knowledge")
        profile_hash = knowledge.get("digest") if isinstance(knowledge, dict) else None
        source_record = {
            "session_id": session_id,
            "session_sha256": _sha256(session_path),
        }
        if isinstance(profile_hash, str) and profile_hash:
            source_record["knowledge_profile_hash"] = profile_hash
        existing_sources = {
            str(item.get("session_sha256", ""))
            for item in staged.get("source_sessions", [])
            if isinstance(item, dict)
        }
        if source_record["session_sha256"] not in existing_sources:
            staged["source_sessions"].append(source_record)
        staged["updated_at"] = datetime.now(timezone.utc).isoformat()
        _write_json_atomic(workspace_path, staged)
        return {
            "ok": True,
            "command": "import-drafts",
            "pack_id": pack_id,
            "imported": len(imported_ids),
            "draft_ids": imported_ids,
        }

    def review_draft(
        self,
        *,
        workspace: Path | str,
        draft_id: str,
        decision: str,
        reviewer: str,
        reason: str,
    ) -> Dict[str, Any]:
        """Record one expert approval or rejection without activating knowledge."""
        if decision not in {"approved", "rejected"}:
            raise KnowledgePackError("review decision must be approved or rejected")
        reviewer = str(reviewer).strip()
        reason = str(reason).strip()
        if not reviewer or not reason:
            raise KnowledgePackError("reviewer and reason are required")
        workspace_path = Path(workspace).resolve()
        staged = _read_json(workspace_path, label="knowledge staging workspace")
        if not isinstance(staged, dict) or staged.get("schema") != "rachel.knowledge-staging.v1":
            raise KnowledgePackError("invalid knowledge staging workspace")
        drafts = staged.get("drafts")
        if not isinstance(drafts, list):
            raise KnowledgePackError("staging workspace has no drafts")
        matches = [
            item
            for item in drafts
            if isinstance(item, dict) and item.get("workspace_draft_id") == draft_id
        ]
        if len(matches) != 1:
            raise KnowledgePackError(f"staging draft not found or ambiguous: {draft_id}")
        draft = matches[0]
        if draft.get("status") != "staging" or "review" in draft:
            raise KnowledgePackError(f"staging draft already reviewed: {draft_id}")
        if draft.get("active") is not False:
            raise KnowledgePackError("staging draft must remain inactive")
        reviewed_at = datetime.now(timezone.utc).isoformat()
        draft["status"] = decision
        draft["active"] = False
        draft["review"] = {
            "decision": decision,
            "reviewer": reviewer,
            "reason": reason,
            "reviewed_at": reviewed_at,
        }
        staged["updated_at"] = reviewed_at
        _write_json_atomic(workspace_path, staged)
        return {
            "ok": True,
            "command": "approve" if decision == "approved" else "reject",
            "draft_id": draft_id,
            "status": decision,
            "active": False,
        }

    def publish(
        self,
        *,
        workspace: Path | str,
        pack_root: Path | str,
        version: str,
        kind: str,
        dependencies: list[str] | None = None,
    ) -> Dict[str, Any]:
        """Publish approved drafts as a new immutable pack version."""
        if not _PACK_PART_RE.fullmatch(str(version)):
            raise KnowledgePackError("knowledge pack version is invalid")
        if kind not in {"team", "project"}:
            raise KnowledgePackError("published pack kind must be team or project")
        workspace_path = Path(workspace).resolve()
        staged = _read_json(workspace_path, label="knowledge staging workspace")
        if not isinstance(staged, dict) or staged.get("schema") != "rachel.knowledge-staging.v1":
            raise KnowledgePackError("invalid knowledge staging workspace")
        pack_id = str(staged.get("pack_id", ""))
        if not _PACK_PART_RE.fullmatch(pack_id):
            raise KnowledgePackError("staging workspace pack id is invalid")

        root = Path(pack_root).resolve()
        final_path = root / pack_id / version
        if final_path.exists():
            raise KnowledgePackError(f"refusing to overwrite published pack: {pack_id}@{version}")
        drafts = staged.get("drafts")
        if not isinstance(drafts, list):
            raise KnowledgePackError("staging workspace has no drafts")
        if any(
            isinstance(draft, dict) and draft.get("status") == "staging"
            for draft in drafts
        ):
            raise KnowledgePackError("all staging drafts must be approved or rejected before publish")
        approved = [
            draft
            for draft in drafts
            if isinstance(draft, dict) and draft.get("status") == "approved"
        ]
        if not approved:
            raise KnowledgePackError("staging workspace has no approved drafts")

        grouped: Dict[str, list[Dict[str, Any]]] = {}
        for draft in approved:
            self._validate_publish_evidence(draft)
            if draft.get("active") is not False:
                raise KnowledgePackError("approved draft must remain inactive before publish")
            if draft.get("target_pack_id") != pack_id:
                raise KnowledgePackError("approved draft target pack mismatch")
            resource = str(draft.get("resource", ""))
            entry = draft.get("entry")
            if not isinstance(entry, dict):
                raise KnowledgePackError("approved draft entry must be an object")
            published_entry = json.loads(json.dumps(entry, ensure_ascii=False))
            if "_knowledge" in published_entry:
                raise KnowledgePackError("draft entry cannot inject publish metadata")
            operation = draft.get("publish") or {"operation": "add"}
            if not isinstance(operation, dict):
                raise KnowledgePackError("draft publish metadata must be an object")
            published_entry["_knowledge"] = json.loads(
                json.dumps(operation, ensure_ascii=False)
            )
            grouped.setdefault(resource, []).append(published_entry)

        dependency_refs = ["rachel.base@1.0.0"]
        for dependency in dependencies or []:
            dependency = str(dependency).strip()
            if dependency and dependency not in dependency_refs:
                dependency_refs.append(dependency)
        root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".rachel_pack_", dir=root) as temp_dir:
            candidate = Path(temp_dir) / "candidate"
            candidate.mkdir()
            manifest_resources: Dict[str, Dict[str, str]] = {}
            for resource in sorted(grouped):
                resource_path = candidate / "resources" / f"{resource}.json"
                _write_json_atomic(
                    resource_path,
                    {"schema": RESOURCE_SCHEMA, "entries": grouped[resource]},
                )
                manifest_resources[resource] = {
                    "path": resource_path.relative_to(candidate).as_posix(),
                    "sha256": _sha256(resource_path),
                }
            manifest_path = candidate / "manifest.json"
            _write_json_atomic(
                manifest_path,
                {
                    "schema": SCHEMA,
                    "id": pack_id,
                    "version": version,
                    "kind": kind,
                    "dependencies": dependency_refs,
                    "resources": manifest_resources,
                },
            )
            self.validate(candidate)
            dependency_packs: list[Path] = []
            for dependency in dependency_refs[1:]:
                dep_id, dep_version = self._split_pack_ref(dependency)
                dep_path = root / dep_id / dep_version
                if not dep_path.is_dir():
                    raise KnowledgePackError(f"publication dependency not found: {dependency}")
                dependency_packs.append(dep_path)
            conflict_result = self.conflicts([*dependency_packs, candidate])
            if not conflict_result["ok"]:
                conflict_types = sorted(
                    {item["type"] for item in conflict_result["conflicts"]}
                )
                raise KnowledgePackError(
                    "publication blocked by knowledge conflicts: "
                    + ", ".join(conflict_types)
                )
            manifest_hash = _sha256(manifest_path)
            final_path.parent.mkdir(parents=True, exist_ok=True)
            candidate.replace(final_path)

        published_at = datetime.now(timezone.utc).isoformat()
        for draft in approved:
            draft["status"] = "published"
            draft["active"] = False
            draft["publication"] = {
                "pack_id": pack_id,
                "version": version,
                "manifest_sha256": manifest_hash,
                "published_at": published_at,
            }
        staged["updated_at"] = published_at
        _write_json_atomic(workspace_path, staged)
        return {
            "ok": True,
            "command": "publish",
            "pack": {"id": pack_id, "version": version, "kind": kind},
            "manifest_sha256": manifest_hash,
            "resource_count": len(grouped),
            "entry_count": len(approved),
        }

    @staticmethod
    def _split_pack_ref(reference: str) -> tuple[str, str]:
        if reference.count("@") != 1:
            raise KnowledgePackError(f"invalid versioned pack dependency: {reference}")
        pack_id, version = reference.split("@", 1)
        if not _PACK_PART_RE.fullmatch(pack_id) or not _PACK_PART_RE.fullmatch(version):
            raise KnowledgePackError(f"invalid versioned pack dependency: {reference}")
        return pack_id, version

    @staticmethod
    def _validate_publish_evidence(draft: Mapping[str, Any]) -> None:
        review = draft.get("review")
        if (
            not isinstance(review, dict)
            or review.get("decision") != "approved"
            or not str(review.get("reviewer", "")).strip()
            or not str(review.get("reason", "")).strip()
        ):
            raise KnowledgePackError("publish requires explicit expert approval")
        source_refs = draft.get("source_refs")
        evidence = draft.get("evidence")
        if (
            not isinstance(source_refs, list)
            or not source_refs
            or not all(
                isinstance(item, dict) and bool(str(item.get("step_id", "")).strip())
                for item in source_refs
            )
        ):
            raise KnowledgePackError("publish requires source step/action references")
        if not isinstance(evidence, list) or not any(
            isinstance(item, dict)
            and item.get("kind") in _ROUTE_EVIDENCE_KINDS
            and any(
                bool(str(item.get(field, "")).strip())
                for field in ("reference", "step_id", "outcome_id")
            )
            for item in evidence
        ):
            raise KnowledgePackError("publish requires route evidence")
        resource = str(draft.get("resource", ""))
        if resource not in _HIGH_EVIDENCE_RESOURCES:
            return
        reaction_instances = [
            item
            for item in evidence
            if isinstance(item, dict) and item.get("kind") == "reaction_instance"
        ]
        if not any(
            KnowledgePackManager._valid_reaction_instance(item)
            for item in reaction_instances
        ):
            raise KnowledgePackError(
                "chemical publication requires a structured reaction instance"
            )
        if not any(
            isinstance(item, dict)
            and item.get("kind") == "validation_result"
            and (
                isinstance(item.get("pass"), bool)
                or bool(str(item.get("status", "")).strip())
            )
            for item in evidence
        ):
            raise KnowledgePackError(
                "chemical publication requires a structured validation result"
            )
        if not any(
            isinstance(item, dict)
            and item.get("kind") in _AUTHORITATIVE_EVIDENCE_KINDS
            and bool(str(item.get("reference", "")).strip())
            for item in evidence
        ):
            raise KnowledgePackError(
                "chemical publication requires experimental, literature, or SOP evidence"
            )

    @staticmethod
    def _valid_reaction_instance(instance: Mapping[str, Any]) -> bool:
        try:
            from rdkit import Chem
        except ImportError as exc:  # pragma: no cover - Rachel runtime requires RDKit
            raise KnowledgePackError("RDKit is required to validate reaction evidence") from exc

        reaction_smiles = instance.get("reaction_smiles")
        if isinstance(reaction_smiles, str) and reaction_smiles.count(">>") == 1:
            left, right = reaction_smiles.split(">>", 1)
            reactants = [item for item in left.split(".") if item]
            products = [item for item in right.split(".") if item]
        else:
            reactants = instance.get("reactants")
            products = instance.get("products")
        if (
            not isinstance(reactants, list)
            or not reactants
            or not isinstance(products, list)
            or not products
            or not all(isinstance(item, str) and item for item in reactants + products)
        ):
            return False
        return all(Chem.MolFromSmiles(item) is not None for item in reactants + products)

    @staticmethod
    def _published_entries(pack: Path | str) -> Dict[str, list[Dict[str, Any]]]:
        manifest_path = _manifest_path(pack)
        manifest = _read_json(manifest_path, label="knowledge manifest")
        result: Dict[str, list[Dict[str, Any]]] = {}
        for resource, spec in manifest["resources"].items():
            overlay = _read_json(
                manifest_path.parent / spec["path"],
                label=f"knowledge resource {resource}",
            )
            result[resource] = list(overlay["entries"])
        return result

    @staticmethod
    def _base_entries(manifest: Mapping[str, Any]) -> Dict[str, Dict[str, Any]]:
        result: Dict[str, Dict[str, Any]] = {}
        for resource, spec in manifest["resources"].items():
            payload = _read_json(
                (_BASE_MANIFEST.parent / spec["path"]).resolve(),
                label=f"base knowledge resource {resource}",
            )
            merge = spec["merge"]
            if merge == "mapping":
                result[resource] = {
                    str(key): value
                    for key, value in payload.items()
                    if not str(key).startswith("__comment")
                }
            elif merge == "catalog":
                entries = payload[spec["entries_path"]]
                id_field = spec.get("id_field", "id")
                result[resource] = {str(item[id_field]): item for item in entries}
            elif merge == "list_catalog":
                id_field = spec.get("id_field", "id")
                result[resource] = {str(item[id_field]): item for item in payload}
            else:  # pragma: no cover - the base manifest is validated at runtime
                raise KnowledgePackError(f"unsupported resource merge mode: {merge}")
        return result

    @staticmethod
    def _entries_by_resource(pack: Path | str) -> Dict[str, Dict[str, str]]:
        result: Dict[str, Dict[str, str]] = {}
        for resource, entries in KnowledgePackManager._published_entries(pack).items():
            result[resource] = {
                str(entry["id"]): json.dumps(
                    entry,
                    ensure_ascii=True,
                    sort_keys=True,
                    separators=(",", ":"),
                )
                for entry in entries
            }
        return result


__all__ = ["KnowledgePackError", "KnowledgePackManager"]
