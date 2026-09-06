"""Audit Rachel terminal leaves against PubChem with conservative chemistry filters.

The script reads a Rachel export directory, terminals.json, or tree.json and
adds a local terminal-material classification before any PubChem lookup. This
prevents salts, counterions, Grignard/organometallic reagents, and inorganic
reagents from being treated as ordinary purchasable building blocks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest"
PUBCHEM_COMPOUND_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
PUBCHEM_VENDOR_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Chemical-Vendors"
PUBCHEM_GHS_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=GHS-Classification"
USER_AGENT = "Rachel-beta-terminal-buyability-audit/0.1"
ALLOWLIST_SOURCE_NAME = "terminal_allowlist.json"
ALLOWLIST_SOURCE_PATH = Path(__file__).resolve().parent.parent / "chem_tools" / ALLOWLIST_SOURCE_NAME

METAL_ATOMIC_NUMBERS = {
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
    30, 37, 38, 39, 40, 44, 45, 46, 47, 48, 50, 55, 56, 74, 75,
    76, 77, 78, 79, 80,
}
ORGANOMETALLIC_METALS = {
    3,   # Li
    11,  # Na
    12,  # Mg
    13,  # Al
    19,  # K
    29,  # Cu
    30,  # Zn
    50,  # Sn
}
HALOGENS = {9, 17, 35, 53}


def canonical_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)


def safe_mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def canonicalize_smiles_text(smiles: str) -> str:
    mol = safe_mol_from_smiles(smiles)
    return canonical_smiles(mol) if mol is not None else ""


@lru_cache(maxsize=1)
def load_terminal_allowlist() -> List[Dict[str, Any]]:
    if not ALLOWLIST_SOURCE_PATH.exists():
        return []
    data = json.loads(ALLOWLIST_SOURCE_PATH.read_text(encoding="utf-8"))
    entries = data.get("entries", [])
    normalized: List[Dict[str, Any]] = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        smiles = str(entry.get("smiles", "")).strip()
        if not smiles:
            continue
        item = dict(entry)
        item["smiles"] = smiles
        item["canonical_smiles"] = canonicalize_smiles_text(smiles)
        item["counts_as_terminal_closure"] = bool(item.get("counts_as_terminal_closure", True))
        item["counts_as_organic_starting_material"] = bool(item.get("counts_as_organic_starting_material", False))
        item["representative_forms"] = [
            dict(form)
            for form in item.get("representative_forms", [])
            if isinstance(form, dict)
        ]
        item["source_reagents"] = [
            dict(reagent)
            for reagent in item.get("source_reagents", [])
            if isinstance(reagent, dict)
        ]
        normalized.append(item)
    return normalized


def find_allowlist_match(smiles: str, local: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    canonical = ""
    if isinstance(local, dict):
        canonical = str(local.get("canonical_smiles", "") or "")
    if not canonical:
        canonical = canonicalize_smiles_text(smiles)
    raw = str(smiles or "").strip()
    query_smiles = str((local or {}).get("query_smiles", "") or "")
    candidates = {raw, canonical, query_smiles}
    for entry in load_terminal_allowlist():
        entry_keys = {entry.get("smiles", ""), entry.get("canonical_smiles", "")}
        if candidates & {key for key in entry_keys if key}:
            return {
                "hit": True,
                "source": ALLOWLIST_SOURCE_NAME,
                "id": str(entry.get("id", "")),
                "label": str(entry.get("label", "")),
                "smiles": str(entry.get("smiles", "")),
                "canonical_smiles": str(entry.get("canonical_smiles", "")),
                "allowlist_class": str(entry.get("allowlist_class", "")),
                "counts_as_terminal_closure": bool(entry.get("counts_as_terminal_closure", True)),
                "counts_as_organic_starting_material": bool(entry.get("counts_as_organic_starting_material", False)),
                "note": str(entry.get("note", "")),
                "evidence_policy": str(entry.get("evidence_policy", "")),
                "representative_forms": list(entry.get("representative_forms", [])),
                "source_reagents": list(entry.get("source_reagents", [])),
                "evidence": {
                    "basis": "no_allowlist_evidence_resolved",
                    "pubchem_cid_closed": False,
                    "vendor_closed": False,
                    "forms": [],
                    "source_reagents": [],
                },
            }
    return {
        "hit": False,
        "source": ALLOWLIST_SOURCE_NAME,
        "id": "",
        "label": "",
        "allowlist_class": "",
        "counts_as_terminal_closure": False,
        "counts_as_organic_starting_material": False,
        "note": "",
        "evidence_policy": "",
        "representative_forms": [],
        "source_reagents": [],
        "evidence": {
            "basis": "no_allowlist_match",
            "pubchem_cid_closed": False,
            "vendor_closed": False,
            "forms": [],
            "source_reagents": [],
        },
    }


def mol_fragments(mol: Chem.Mol) -> List[Chem.Mol]:
    fragments = []
    for frag in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False):
        try:
            Chem.SanitizeMol(frag)
        except Exception:
            pass
        fragments.append(frag)
    return fragments


def atom_counts(mol: Chem.Mol) -> Dict[str, int]:
    heavy = 0
    carbons = 0
    metals = 0
    halogens = 0
    charge = 0
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num > 1:
            heavy += 1
        if atomic_num == 6:
            carbons += 1
        if atomic_num in METAL_ATOMIC_NUMBERS:
            metals += 1
        if atomic_num in HALOGENS:
            halogens += 1
        charge += atom.GetFormalCharge()
    return {
        "heavy_atoms": heavy,
        "carbons": carbons,
        "metals": metals,
        "halogens": halogens,
        "formal_charge": charge,
    }


def has_direct_carbon_metal_bond(mol: Chem.Mol) -> bool:
    for bond in mol.GetBonds():
        nums = {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
        if 6 in nums and nums.intersection(ORGANOMETALLIC_METALS):
            return True
    return False


def has_carbanion(mol: Chem.Mol) -> bool:
    return any(atom.GetAtomicNum() == 6 and atom.GetFormalCharge() < 0 for atom in mol.GetAtoms())


def uncharge_fragment(mol: Chem.Mol) -> Chem.Mol:
    try:
        uncharger = rdMolStandardize.Uncharger()
        uncharged = uncharger.uncharge(Chem.Mol(mol))
        Chem.SanitizeMol(uncharged)
        return uncharged
    except Exception:
        return mol


def choose_largest_organic_fragment(fragments: Iterable[Chem.Mol]) -> Optional[Chem.Mol]:
    organic = [frag for frag in fragments if atom_counts(frag)["carbons"] > 0]
    if not organic:
        return None
    return max(
        organic,
        key=lambda frag: (
            atom_counts(frag)["carbons"],
            atom_counts(frag)["heavy_atoms"],
            len(canonical_smiles(frag)),
        ),
    )


def classify_terminal_smiles(smiles: str) -> Dict[str, Any]:
    """Return a conservative local classification and PubChem query policy."""
    raw_smiles = str(smiles or "").strip()
    result: Dict[str, Any] = {
        "input_smiles": raw_smiles,
        "valid_smiles": False,
        "canonical_smiles": "",
        "query_smiles": "",
        "query_mode": "skip_invalid",
        "local_classification": "invalid_smiles",
        "terminal_material_policy": "manual_review",
        "risk_flags": [],
        "normalization_notes": [],
        "fragments": [],
    }
    mol = safe_mol_from_smiles(raw_smiles)
    if mol is None:
        result["risk_flags"].append("invalid_smiles")
        return result

    result["valid_smiles"] = True
    result["canonical_smiles"] = canonical_smiles(mol)
    counts = atom_counts(mol)
    fragments = mol_fragments(mol)
    result["fragments"] = [
        {
            "smiles": canonical_smiles(frag),
            **atom_counts(frag),
        }
        for frag in fragments
    ]

    has_metal = counts["metals"] > 0
    organic_parent = choose_largest_organic_fragment(fragments)
    organic_parent_counts = atom_counts(organic_parent) if organic_parent is not None else {}
    direct_organometallic = has_direct_carbon_metal_bond(mol)
    metal_same_fragment_as_carbon = bool(
        organic_parent is not None and organic_parent_counts.get("metals", 0) > 0
    )
    ionic_organometallic = bool(has_metal and organic_parent is not None and has_carbanion(organic_parent))

    if direct_organometallic or metal_same_fragment_as_carbon or ionic_organometallic:
        result.update({
            "query_mode": "skip_reactive_reagent",
            "local_classification": "reactive_organometallic_reagent",
            "terminal_material_policy": "not_suitable_as_terminal_material",
        })
        if direct_organometallic:
            result["risk_flags"].append("direct_carbon_metal_bond")
        if metal_same_fragment_as_carbon:
            result["risk_flags"].append("metal_in_organic_fragment")
        if ionic_organometallic:
            result["risk_flags"].append("carbanion_with_metal_counterion")
        result["normalization_notes"].append(
            "Organolithium, Grignard, organozinc, organoaluminum, organocopper, or organotin-like species should be treated as reagents or generated in situ, not as ordinary terminal building blocks."
        )
        return result

    if counts["carbons"] == 0:
        result.update({
            "query_mode": "skip_inorganic_or_reagent",
            "local_classification": "inorganic_or_elemental_reagent",
            "terminal_material_policy": "reagent_not_route_terminal",
        })
        result["risk_flags"].append("no_carbon_atoms")
        result["normalization_notes"].append(
            "Elemental/inorganic reagents may be purchasable but should not close an organic precursor tree as building-block terminals."
        )
        return result

    if has_metal and organic_parent is not None:
        parent = uncharge_fragment(organic_parent)
        parent_smiles = canonical_smiles(parent)
        result.update({
            "query_smiles": parent_smiles,
            "query_mode": "query_parent_from_salt",
            "local_classification": "metal_counterion_or_salt",
            "terminal_material_policy": "parent_query_manual_salt_review",
        })
        result["risk_flags"].append("metal_counterion_stripped")
        if parent_smiles != canonical_smiles(organic_parent):
            result["risk_flags"].append("parent_neutralized")
        result["normalization_notes"].append(
            "Metal counterion stripped before PubChem parent query; manually decide whether salt form or neutral parent is the real purchasable terminal."
        )
        return result

    if len(fragments) > 1 and organic_parent is not None:
        parent = uncharge_fragment(organic_parent)
        parent_smiles = canonical_smiles(parent)
        result.update({
            "query_smiles": parent_smiles,
            "query_mode": "query_largest_organic_fragment",
            "local_classification": "organic_mixture_or_salt",
            "terminal_material_policy": "parent_query_manual_mixture_review",
        })
        result["risk_flags"].append("multiple_fragments")
        result["normalization_notes"].append(
            "Multiple fragments detected; PubChem query uses the largest organic parent rather than the full mixture."
        )
        return result

    parent = uncharge_fragment(mol)
    parent_smiles = canonical_smiles(parent)
    result.update({
        "query_smiles": parent_smiles,
        "query_mode": "query_direct",
        "local_classification": "ordinary_organic_candidate",
        "terminal_material_policy": "direct_terminal_candidate",
    })
    if parent_smiles != result["canonical_smiles"]:
        result["risk_flags"].append("neutralized_for_query")
        result["normalization_notes"].append("Formal charge neutralized before PubChem query.")
    return result


class PubChemClient:
    def __init__(
        self,
        *,
        cache_dir: Optional[Path] = None,
        timeout: float = 20.0,
        pause_seconds: float = 0.2,
        offline: bool = False,
    ) -> None:
        self.cache_dir = cache_dir
        self.timeout = timeout
        self.pause_seconds = pause_seconds
        self.offline = offline
        if cache_dir:
            cache_dir.mkdir(parents=True, exist_ok=True)

    def _cache_path(self, key: str) -> Optional[Path]:
        if self.cache_dir is None:
            return None
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.json"

    def request_json(self, key: str, url: str, *, data: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
        cache_path = self._cache_path(key)
        if cache_path and cache_path.exists():
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if not cached.get("error"):
                request_meta = dict(cached.get("_rachel_request") or {})
                request_meta["from_cache"] = True
                cached["_rachel_request"] = request_meta
                return cached
        if self.offline:
            return {"offline": True}

        encoded_data = urllib.parse.urlencode(data).encode("utf-8") if data else None
        request = urllib.request.Request(
            url,
            data=encoded_data,
            headers={"User-Agent": USER_AGENT},
            method="POST" if data else "GET",
        )
        payload: Dict[str, Any] = {}
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    payload = json.load(response)
                break
            except urllib.error.HTTPError as exc:
                if exc.code == 404:
                    payload = {"not_found": True, "status": exc.code, "reason": str(exc)}
                    break
                payload = {"error": "http_error", "status": exc.code, "reason": str(exc)}
            except Exception as exc:
                payload = {"error": type(exc).__name__, "reason": str(exc)}
            time.sleep(min(1.0, 0.2 * (attempt + 1)))

        payload["_rachel_request"] = {
            "url": url,
            "retrieved_at": datetime.now(timezone.utc).isoformat(),
            "from_cache": False,
        }

        if cache_path and not payload.get("error"):
            cache_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        if self.pause_seconds:
            time.sleep(self.pause_seconds)
        return payload

    def cid_lookup(self, smiles: str) -> Dict[str, Any]:
        url = f"{PUBCHEM_BASE}/pug/compound/smiles/cids/JSON"
        payload = self.request_json(
            f"smiles_to_cids:{smiles}",
            url,
            data={"smiles": smiles},
        )
        if payload.get("offline"):
            return {"status": "not_queried", "cids": [], "error": ""}
        if payload.get("error"):
            return {
                "status": "lookup_error",
                "cids": [],
                "error": payload.get("reason") or payload.get("error", ""),
            }
        if payload.get("not_found") or payload.get("Fault"):
            return {"status": "no_hit", "cids": [], "error": ""}
        cids = payload.get("IdentifierList", {}).get("CID", [])
        normalized = [
            int(cid)
            for cid in cids
            if (isinstance(cid, int) or str(cid).isdigit()) and int(cid) > 0
        ]
        return {"status": "hit" if normalized else "no_hit", "cids": normalized, "error": ""}

    def name_lookup(self, name: str) -> Dict[str, Any]:
        encoded = urllib.parse.quote(name.strip())
        url = f"{PUBCHEM_BASE}/pug/compound/name/{encoded}/cids/JSON"
        payload = self.request_json(f"name_to_cids:{name}", url)
        if payload.get("offline"):
            return {"status": "not_queried", "cids": [], "error": ""}
        if payload.get("error"):
            return {
                "status": "lookup_error",
                "cids": [],
                "error": payload.get("reason") or payload.get("error", ""),
            }
        if payload.get("not_found") or payload.get("Fault"):
            return {"status": "no_hit", "cids": [], "error": ""}
        cids = payload.get("IdentifierList", {}).get("CID", [])
        normalized = [
            int(cid)
            for cid in cids
            if (isinstance(cid, int) or str(cid).isdigit()) and int(cid) > 0
        ]
        return {"status": "hit" if normalized else "no_hit", "cids": normalized, "error": ""}

    def cids_from_smiles(self, smiles: str) -> List[int]:
        return list(self.cid_lookup(smiles).get("cids", []))

    def properties(self, cid: int) -> Dict[str, Any]:
        props = (
            "MolecularFormula,Charge,CanonicalSMILES,ConnectivitySMILES,"
            "IsomericSMILES,InChIKey,IUPACName"
        )
        url = f"{PUBCHEM_BASE}/pug/compound/cid/{cid}/property/{props}/JSON"
        payload = self.request_json(f"cid_properties:{cid}", url)
        rows = payload.get("PropertyTable", {}).get("Properties", [])
        return rows[0] if rows else {}

    def vendor_summary(self, cid: int) -> Dict[str, Any]:
        data_url = f"{PUBCHEM_BASE}/pug_view/data/compound/{cid}/JSON?heading=Chemical%20Vendors"
        category_url = f"{PUBCHEM_BASE}/pug_view/categories/compound/{cid}/JSON?heading=Chemical%20Vendors"
        data_payload = self.request_json(f"cid_vendor_flag:{cid}", data_url)
        category_payload = self.request_json(f"cid_vendor_categories:{cid}", category_url)
        if data_payload.get("offline") or category_payload.get("offline"):
            return {
                "queried": False,
                "vendor_flag": None,
                "vendor_count": None,
                "vendor_names": [],
                "vendor_rows": [],
                "vendor_catalog_status": "not_queried",
                "vendor_url": "",
            }
        rows = _parse_vendor_rows(category_payload)
        all_names = sorted({row["source_name"] for row in rows if row.get("source_name")})
        category_error = category_payload.get("reason") or category_payload.get("error", "")
        aggregate_flag = _contains_true_boolean(data_payload)
        if rows:
            catalog_status = "hit"
        elif category_error:
            catalog_status = "lookup_error"
        elif aggregate_flag:
            catalog_status = "aggregate_flag_only"
        else:
            catalog_status = "no_hit"
        request_meta = category_payload.get("_rachel_request") or {}
        return {
            "queried": True,
            "vendor_flag": bool(rows),
            "vendor_count": len(all_names),
            "vendor_names": all_names[:25],
            "vendor_listing_count": len(rows),
            "vendor_rows": rows[:50],
            "vendor_rows_truncated": len(rows) > 50,
            "vendor_catalog_status": catalog_status,
            "aggregate_vendor_flag": bool(aggregate_flag),
            "vendor_url": PUBCHEM_VENDOR_URL.format(cid=cid),
            "retrieved_at": request_meta.get("retrieved_at", ""),
            "from_cache": bool(request_meta.get("from_cache")),
            "error": str(category_error or ""),
            "note": (
                "PubChem Chemical Vendors rows are depositor catalog evidence, "
                "not live inventory, price, lead time, or procurement confirmation."
            ),
        }

    def ghs_summary(self, cid: int) -> Dict[str, Any]:
        """Return source-attributed PubChem GHS classifications for one CID."""
        url = f"{PUBCHEM_BASE}/pug_view/data/compound/{cid}/JSON?heading=GHS%20Classification"
        payload = self.request_json(f"cid_ghs:{cid}", url)
        result: Dict[str, Any] = {
            "queried": not self.offline,
            "state": "not_queried" if self.offline else "no_classification_returned",
            "cid": cid,
            "url": PUBCHEM_GHS_URL.format(cid=cid),
            "signal_words": [],
            "pictograms": [],
            "hazard_statements": [],
            "precautionary_statement_codes": [],
            "classifications": [],
            "note": (
                "PubChem aggregates classifications from multiple sources. "
                "These records are not a substitute for the product- and supplier-specific SDS."
            ),
        }
        request_meta = payload.get("_rachel_request") or {}
        result["retrieved_at"] = request_meta.get("retrieved_at", "")
        result["from_cache"] = bool(request_meta.get("from_cache"))
        if payload.get("offline"):
            return result
        if payload.get("error"):
            result["state"] = "lookup_error"
            result["error"] = payload.get("reason") or payload.get("error", "")
            return result
        if payload.get("not_found") or payload.get("Fault"):
            return result

        classifications = _parse_ghs_classifications(payload)
        result["classifications"] = classifications
        result["signal_words"] = _unique_strings(
            value for item in classifications for value in item.get("signal_words", [])
        )
        result["pictograms"] = _unique_strings(
            value for item in classifications for value in item.get("pictograms", [])
        )
        result["hazard_statements"] = _unique_strings(
            value
            for item in classifications
            for value in item.get("hazard_statements", [])
        )
        result["precautionary_statement_codes"] = _unique_strings(
            value
            for item in classifications
            for value in item.get("precautionary_statement_codes", [])
        )
        if classifications:
            result["state"] = "classification_data_found"
        return result


def _contains_true_boolean(value: Any) -> bool:
    if isinstance(value, dict):
        return any(_contains_true_boolean(item) for item in value.values())
    if isinstance(value, list):
        return any(_contains_true_boolean(item) for item in value)
    return value is True


def _collect_vendor_like_names(value: Any) -> List[str]:
    names: List[str] = []
    if isinstance(value, dict):
        for key, item in value.items():
            if key in {"SourceName", "Name", "RegistryName"} and isinstance(item, str):
                if item and item.lower() not in {"pubchem", "compound"}:
                    names.append(item)
            names.extend(_collect_vendor_like_names(item))
    elif isinstance(value, list):
        for item in value:
            names.extend(_collect_vendor_like_names(item))
    return names


def _parse_vendor_rows(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    source_categories = payload.get("SourceCategories") or {}
    categories = source_categories.get("Categories") or []
    if isinstance(categories, dict):
        categories = [categories]
    rows: List[Dict[str, Any]] = []
    for category in categories:
        if not isinstance(category, dict) or category.get("Category") != "Chemical Vendors":
            continue
        for source in category.get("Sources", []) or []:
            if not isinstance(source, dict):
                continue
            sid = _positive_int(source.get("SID"))
            rows.append({
                "sid": sid,
                "source_name": str(source.get("SourceName", "") or "").strip(),
                "source_url": str(source.get("SourceURL", "") or "").strip(),
                "source_detail_url": str(source.get("SourceDetail", "") or "").strip(),
                "registry_id": str(source.get("RegistryID", "") or "").strip(),
                "source_record_url": str(source.get("SourceRecordURL", "") or "").strip(),
                "source_categories": [
                    str(value)
                    for value in source.get("SourceCategories", []) or []
                    if str(value).strip()
                ],
            })
    rows.sort(key=lambda item: (
        item.get("source_name", "").lower(),
        item.get("registry_id", "").lower(),
        int(item.get("sid") or 0),
    ))
    return rows


def _unique_strings(values: Iterable[Any]) -> List[str]:
    unique: List[str] = []
    seen = set()
    for value in values:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            unique.append(text)
    return unique


def _find_sections(value: Any, heading: str) -> List[Dict[str, Any]]:
    sections: List[Dict[str, Any]] = []
    if isinstance(value, dict):
        if value.get("TOCHeading") == heading:
            sections.append(value)
        for item in value.get("Section", []) or []:
            sections.extend(_find_sections(item, heading))
    elif isinstance(value, list):
        for item in value:
            sections.extend(_find_sections(item, heading))
    return sections


def _information_strings(information: Dict[str, Any]) -> List[str]:
    value = information.get("Value") or {}
    strings: List[str] = []
    for item in value.get("StringWithMarkup", []) or []:
        if isinstance(item, dict):
            strings.append(item.get("String", ""))
    if isinstance(value.get("String"), str):
        strings.append(value["String"])
    return _unique_strings(strings)


def _pictogram_names(information: Dict[str, Any]) -> List[str]:
    value = information.get("Value") or {}
    names: List[str] = []
    for item in value.get("StringWithMarkup", []) or []:
        if not isinstance(item, dict):
            continue
        for markup in item.get("Markup", []) or []:
            if not isinstance(markup, dict) or markup.get("Type") != "Icon":
                continue
            names.append(markup.get("Extra", ""))
    return _unique_strings(names)


def _parse_ghs_classifications(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    record = payload.get("Record") or {}
    references = {
        int(item.get("ReferenceNumber")): item
        for item in record.get("Reference", []) or []
        if isinstance(item, dict) and _positive_int(item.get("ReferenceNumber"))
    }
    grouped: Dict[int, Dict[str, Any]] = {}
    for section in _find_sections(record, "GHS Classification"):
        for information in section.get("Information", []) or []:
            if not isinstance(information, dict):
                continue
            reference_number = _positive_int(information.get("ReferenceNumber")) or 0
            group = grouped.setdefault(
                reference_number,
                {
                    "reference_number": reference_number or None,
                    "source_name": "",
                    "source_id": "",
                    "source_url": "",
                    "license_url": "",
                    "signal_words": [],
                    "pictograms": [],
                    "hazard_statements": [],
                    "precautionary_statement_codes": [],
                },
            )
            name = str(information.get("Name", "") or "")
            if name == "Pictogram(s)":
                group["pictograms"].extend(_pictogram_names(information))
            elif name == "Signal":
                group["signal_words"].extend(_information_strings(information))
            elif name == "GHS Hazard Statements":
                group["hazard_statements"].extend(_information_strings(information))
            elif name == "Precautionary Statement Codes":
                group["precautionary_statement_codes"].extend(
                    _information_strings(information)
                )

    classifications: List[Dict[str, Any]] = []
    for reference_number, group in grouped.items():
        reference = references.get(reference_number, {})
        group["source_name"] = str(reference.get("SourceName", "") or "")
        group["source_id"] = str(reference.get("SourceID", "") or "")
        group["source_url"] = str(reference.get("URL", "") or "")
        group["license_note"] = str(reference.get("LicenseNote", "") or "")
        group["license_url"] = str(reference.get("LicenseURL", "") or "")
        for key in (
            "signal_words",
            "pictograms",
            "hazard_statements",
            "precautionary_statement_codes",
        ):
            group[key] = _unique_strings(group[key])
        if group["signal_words"] or group["pictograms"] or group["hazard_statements"]:
            classifications.append(group)
    classifications.sort(key=lambda item: int(item.get("reference_number") or 0))
    return classifications


def _positive_int(value: Any) -> Optional[int]:
    if isinstance(value, int) and value > 0:
        return value
    if isinstance(value, str) and value.isdigit() and int(value) > 0:
        return int(value)
    return None


def lookup_evidence_material(
    material: Dict[str, Any],
    client: PubChemClient,
    *,
    include_vendors: bool,
) -> Dict[str, Any]:
    """Resolve an allowlist real-world material to PubChem CID/vendor evidence."""
    name = str(material.get("name", "") or material.get("label", "") or "").strip()
    smiles = str(material.get("smiles", "") or "").strip()
    role = str(material.get("role", "") or "").strip()
    result: Dict[str, Any] = {
        "name": name,
        "smiles": smiles,
        "role": role,
        "note": str(material.get("note", "") or ""),
        "query_type": "",
        "query_value": "",
        "cid_lookup": {"status": "not_queried", "cids": [], "error": ""},
        "cids": [],
        "cid": None,
        "cid_url": "",
        "cid_hit": False,
        "vendor_evidence": {
            "queried": False,
            "vendor_flag": None,
            "vendor_count": None,
            "vendor_names": [],
            "vendor_url": "",
        },
        "vendor_hit": False,
        "properties": {},
    }
    explicit_cid = _positive_int(material.get("cid"))
    if client.offline:
        return result
    if explicit_cid is not None:
        result["query_type"] = "cid"
        result["query_value"] = str(explicit_cid)
        result["cid_lookup"] = {"status": "hit", "cids": [explicit_cid], "error": ""}
        cids = [explicit_cid]
    elif smiles:
        result["query_type"] = "smiles"
        result["query_value"] = smiles
        lookup = client.cid_lookup(smiles)
        result["cid_lookup"] = lookup
        cids = list(lookup.get("cids", []))
    elif name:
        result["query_type"] = "name"
        result["query_value"] = name
        lookup = client.name_lookup(name)
        result["cid_lookup"] = lookup
        cids = list(lookup.get("cids", []))
    else:
        cids = []

    positive_cids = [
        int(cid)
        for cid in cids
        if (isinstance(cid, int) or str(cid).isdigit()) and int(cid) > 0
    ]
    result["cids"] = positive_cids[:10]
    result["cid_hit"] = bool(positive_cids)
    if positive_cids:
        best_cid = positive_cids[0]
        result["cid"] = best_cid
        result["cid_url"] = PUBCHEM_COMPOUND_URL.format(cid=best_cid)
        result["properties"] = client.properties(best_cid)
        if include_vendors:
            vendor = client.vendor_summary(best_cid)
            result["vendor_evidence"] = vendor
            result["vendor_hit"] = bool(vendor.get("vendor_flag") is True)
    return result


def resolve_allowlist_evidence(
    allowlist: Dict[str, Any],
    client: PubChemClient,
    *,
    include_vendors: bool,
) -> Dict[str, Any]:
    """Attach CID/vendor evidence for an allowlisted reagent or placeholder leaf."""
    if not allowlist.get("hit"):
        return allowlist
    resolved = dict(allowlist)
    forms = [
        lookup_evidence_material(form, client, include_vendors=include_vendors)
        for form in allowlist.get("representative_forms", [])
        if isinstance(form, dict)
    ]
    sources = [
        lookup_evidence_material(reagent, client, include_vendors=include_vendors)
        for reagent in allowlist.get("source_reagents", [])
        if isinstance(reagent, dict)
    ]

    form_cid_closed = any(form.get("cid_hit") for form in forms)
    form_vendor_closed = any(form.get("vendor_hit") for form in forms)
    source_cid_closed = bool(sources) and all(source.get("cid_hit") for source in sources)
    source_vendor_closed = bool(sources) and all(source.get("vendor_hit") for source in sources)

    if form_vendor_closed:
        basis = "allowlist_representative_form"
        cid_closed = True
        vendor_closed = True
    elif source_vendor_closed:
        basis = "allowlist_in_situ_sources"
        cid_closed = True
        vendor_closed = True
    elif form_cid_closed:
        basis = "allowlist_representative_form_cid_only"
        cid_closed = True
        vendor_closed = False
    elif source_cid_closed:
        basis = "allowlist_in_situ_sources_cid_only"
        cid_closed = True
        vendor_closed = False
    else:
        basis = "allowlist_without_pubchem_evidence"
        cid_closed = False
        vendor_closed = False

    resolved["evidence"] = {
        "basis": basis,
        "pubchem_cid_closed": cid_closed,
        "vendor_closed": vendor_closed,
        "forms": forms,
        "source_reagents": sources,
        "note": "Allowlist evidence maps an abstract reagent/ion/placeholder leaf to real-world purchasable or in-situ source materials.",
    }
    return resolved


def load_terminal_records(path: Path) -> Tuple[List[Dict[str, Any]], Path]:
    resolved = path.resolve()
    if resolved.is_dir():
        candidates = [resolved / "terminals.json", resolved / "export" / "terminals.json"]
        terminal_path = next((item for item in candidates if item.exists()), None)
        if terminal_path is None:
            raise FileNotFoundError(f"No terminals.json found under {resolved}")
    else:
        terminal_path = resolved

    data = json.loads(terminal_path.read_text(encoding="utf-8"))
    if isinstance(data, list):
        return [dict(item) for item in data if isinstance(item, dict)], terminal_path
    if isinstance(data, dict) and isinstance(data.get("terminals"), list):
        return [dict(item) for item in data["terminals"] if isinstance(item, dict)], terminal_path
    if isinstance(data, dict) and "molecule_nodes" in data:
        records = []
        for node_id, node in data.get("molecule_nodes", {}).items():
            if not isinstance(node, dict):
                continue
            if node.get("role") == "terminal" or (node.get("complexity") or {}).get("is_terminal"):
                records.append({
                    "node_id": node.get("node_id", node_id),
                    "smiles": node.get("smiles", ""),
                    "cs_score": (node.get("complexity") or {}).get("cs_score"),
                    "classification": (node.get("complexity") or {}).get("classification", ""),
                })
        return records, terminal_path
    raise ValueError(f"Unsupported terminal input JSON structure: {terminal_path}")


def audit_record(
    record: Dict[str, Any],
    client: PubChemClient,
    *,
    include_vendors: bool,
    query_reagents: bool,
) -> Dict[str, Any]:
    smiles = str(record.get("smiles", "")).strip()
    local = classify_terminal_smiles(smiles)
    allowlist = find_allowlist_match(smiles, local)
    allowlist = resolve_allowlist_evidence(
        allowlist,
        client,
        include_vendors=include_vendors,
    )
    pubchem: Dict[str, Any] = {
        "queried": False,
        "query_smiles": local.get("query_smiles", ""),
        "cid_lookup": {"status": "not_queried", "cids": [], "error": ""},
        "cids": [],
        "best_cid": None,
        "best_cid_url": "",
        "properties": {},
        "vendor_evidence": {
            "queried": False,
            "vendor_flag": None,
            "vendor_count": None,
            "vendor_names": [],
            "vendor_url": "",
        },
    }

    can_query = (
        not client.offline
        and bool(local.get("query_smiles"))
        and local.get("query_mode", "").startswith("query_")
    )
    if query_reagents and local.get("canonical_smiles") and local.get("query_mode", "").startswith("skip_"):
        pubchem["query_smiles"] = local["canonical_smiles"]
        can_query = not client.offline

    if can_query:
        pubchem["queried"] = True
        lookup = client.cid_lookup(str(pubchem["query_smiles"]))
        pubchem["cid_lookup"] = lookup
        cids = list(lookup.get("cids", []))
        pubchem["cids"] = cids[:10]
        if cids:
            best_cid = cids[0]
            pubchem["best_cid"] = best_cid
            pubchem["best_cid_url"] = PUBCHEM_COMPOUND_URL.format(cid=best_cid)
            pubchem["properties"] = client.properties(best_cid)
            if include_vendors:
                pubchem["vendor_evidence"] = client.vendor_summary(best_cid)

    metrics = build_pubchem_metrics(local, pubchem, allowlist)
    decision = decide_buyability(local, pubchem, metrics, allowlist)
    return {
        "node_id": record.get("node_id", ""),
        "smiles": smiles,
        "cs_score": record.get("cs_score"),
        "rachel_classification": record.get("classification", ""),
        "local": local,
        "allowlist": allowlist,
        "pubchem": pubchem,
        "pubchem_metrics": metrics,
        "buyability_decision": decision,
    }


def build_pubchem_metrics(
    local: Dict[str, Any],
    pubchem: Dict[str, Any],
    allowlist: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Project PubChem lookup into two explicit evidence indicators.

    cid_hit means PubChem has a public compound record for the query. vendor_hit
    is a stronger Chemical Vendors signal, but it is not live inventory, price,
    or shipping evidence.
    """
    cids = [int(cid) for cid in pubchem.get("cids", []) if isinstance(cid, int) or str(cid).isdigit()]
    vendor = pubchem.get("vendor_evidence") or {}
    vendor_flag = vendor.get("vendor_flag")
    queried = bool(pubchem.get("queried"))
    vendor_queried = bool(vendor.get("queried"))
    lookup = pubchem.get("cid_lookup") or {}
    lookup_status = str(lookup.get("status", "") or "")
    cid_hit = bool(cids)
    vendor_hit = bool(vendor_flag is True)
    allowlist = allowlist or {}
    allowlist_hit = bool(allowlist.get("hit"))
    allowlist_evidence = allowlist.get("evidence") if isinstance(allowlist.get("evidence"), dict) else {}
    evidence_cid_closed = bool(allowlist_hit and allowlist_evidence.get("pubchem_cid_closed"))
    evidence_vendor_closed = bool(allowlist_hit and allowlist_evidence.get("vendor_closed"))
    allowlist_basis = str(allowlist_evidence.get("basis", "") or "")
    pubchem_cid_closed = bool(cid_hit or evidence_cid_closed)
    vendor_closed = bool(vendor_hit or evidence_vendor_closed)
    if vendor_hit:
        closure_basis = "chemical_vendor"
    elif evidence_vendor_closed:
        closure_basis = allowlist_basis or "allowlist_vendor_evidence"
    elif cid_hit:
        closure_basis = "pubchem_public_record_only"
    elif evidence_cid_closed:
        closure_basis = allowlist_basis or "allowlist_pubchem_cid_evidence"
        if not closure_basis.endswith("_cid_only"):
            closure_basis = f"{closure_basis}_cid_only"
    elif lookup_status == "lookup_error":
        closure_basis = "pubchem_lookup_error"
    else:
        closure_basis = "not_closed"
    if lookup_status == "lookup_error":
        cid_status = "lookup_error"
    elif not queried:
        cid_status = "not_queried"
    elif cid_hit:
        cid_status = "hit"
    else:
        cid_status = "no_hit"
    if not vendor_queried:
        vendor_status = "not_queried"
    elif vendor_hit:
        vendor_status = "hit"
    else:
        vendor_status = "no_hit"

    return {
        "query_mode": local.get("query_mode", ""),
        "query_smiles": pubchem.get("query_smiles", ""),
        "cid_hit": cid_hit,
        "public_record": cid_hit,
        "cid_status": cid_status,
        "cid_lookup_error": lookup.get("error", ""),
        "cid_count": len(cids),
        "best_cid": pubchem.get("best_cid"),
        "vendor_hit": vendor_hit,
        "vendor_status": vendor_status,
        "vendor_count": vendor.get("vendor_count"),
        "vendor_names": list(vendor.get("vendor_names") or []),
        "best_cid_url": PUBCHEM_COMPOUND_URL.format(cid=pubchem.get("best_cid")) if pubchem.get("best_cid") else "",
        "cid_urls": [
            PUBCHEM_COMPOUND_URL.format(cid=cid)
            for cid in cids[:10]
        ],
        "vendor_url": vendor.get("vendor_url", ""),
        "counts_as_organic_starting_material": bool(
            allowlist.get("counts_as_organic_starting_material", False)
        ) if allowlist_hit else bool(vendor_hit),
        "pubchem_cid_closed": pubchem_cid_closed,
        "vendor_closed": vendor_closed,
        "closure_basis": closure_basis,
        "closure_evidence": {
            "basis": closure_basis,
            "direct_query": {
                "cid_hit": cid_hit,
                "vendor_hit": vendor_hit,
                "query_smiles": pubchem.get("query_smiles", ""),
                "cids": cids[:10],
                "cid_urls": [
                    PUBCHEM_COMPOUND_URL.format(cid=cid)
                    for cid in cids[:10]
                ],
                "vendor_url": vendor.get("vendor_url", ""),
            },
            "allowlist": allowlist_evidence if allowlist_hit else {},
        },
    }


def decide_buyability(
    local: Dict[str, Any],
    pubchem: Dict[str, Any],
    metrics: Dict[str, Any],
    allowlist: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    policy = local.get("terminal_material_policy", "")
    query_mode = local.get("query_mode", "")
    allowlist = allowlist or {}
    closure_evidence = metrics.get("closure_evidence", {})
    allowlist_evidence = closure_evidence.get("allowlist", {}) if isinstance(closure_evidence, dict) else {}
    evidence_vendor_closed = bool(allowlist_evidence.get("vendor_closed"))
    evidence_cid_closed = bool(allowlist_evidence.get("pubchem_cid_closed"))
    if evidence_vendor_closed and policy in {
        "not_suitable_as_terminal_material",
        "reagent_not_route_terminal",
    }:
        return {
            "state": "vendor_evidence_allowlist_reagent_form",
            "confidence": "high",
            "reason": "The abstract reagent/ion/metal leaf is mapped to allowlisted real-world form or source reagents with PubChem Chemical Vendors evidence.",
            "terminal_material_policy": policy,
            "allowlist_id": allowlist.get("id", ""),
            "allowlist_class": allowlist.get("allowlist_class", ""),
            "closure_basis": metrics.get("closure_basis", ""),
        }
    if evidence_cid_closed and policy in {
        "not_suitable_as_terminal_material",
        "reagent_not_route_terminal",
    }:
        return {
            "state": "cid_evidence_allowlist_reagent_form_manual_vendor_check",
            "confidence": "medium",
            "reason": "The abstract reagent/ion/metal leaf maps to PubChem CID evidence, but Chemical Vendors evidence was not found for the mapped form/source set.",
            "terminal_material_policy": policy,
            "allowlist_id": allowlist.get("id", ""),
            "allowlist_class": allowlist.get("allowlist_class", ""),
            "closure_basis": metrics.get("closure_basis", ""),
        }
    if policy in {"not_suitable_as_terminal_material", "reagent_not_route_terminal"}:
        return {
            "state": "not_route_terminal_material",
            "confidence": "high",
            "reason": "Local chemistry filter classifies this leaf as a reagent/counter-reagent rather than an ordinary purchasable route terminal.",
            "terminal_material_policy": policy,
        }
    if not pubchem.get("queried"):
        return {
            "state": "not_queried_manual_review",
            "confidence": "medium",
            "reason": "No PubChem query was run for this local classification.",
        }
    if metrics.get("cid_status") == "lookup_error":
        return {
            "state": "pubchem_lookup_error_manual_retry",
            "confidence": "unknown",
            "reason": "PubChem CID lookup failed; retry before treating this as no public record.",
        }
    if not metrics.get("cid_hit"):
        return {
            "state": "no_public_pubchem_record",
            "confidence": "medium",
            "reason": "PubChem did not return a CID for the query SMILES.",
        }

    if metrics.get("vendor_hit") is True and query_mode == "query_direct":
        return {
            "state": "vendor_evidence_direct",
            "confidence": "medium_high",
            "reason": "PubChem returned a CID and Chemical Vendors evidence for the direct terminal query.",
        }
    if metrics.get("vendor_hit") is True:
        return {
            "state": "vendor_evidence_parent_manual_form_check",
            "confidence": "medium",
            "reason": "PubChem vendor evidence exists for a normalized parent; confirm salt/mixture/form is acceptable.",
        }
    return {
        "state": "public_record_only_manual_vendor_check",
        "confidence": "medium_low",
        "reason": "PubChem returned a CID, but no Chemical Vendors evidence was parsed in this audit.",
    }


def summarize(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    by_state: Dict[str, int] = {}
    by_local: Dict[str, int] = {}
    by_cid_status: Dict[str, int] = {}
    by_vendor_status: Dict[str, int] = {}
    pubchem_cid_closed = 0
    vendor_closed = 0
    organic_starting_material_closures = 0
    for item in results:
        state = item.get("buyability_decision", {}).get("state", "unknown")
        local = item.get("local", {}).get("local_classification", "unknown")
        metrics = item.get("pubchem_metrics", {})
        cid_status = metrics.get("cid_status", "unknown")
        vendor_status = metrics.get("vendor_status", "unknown")
        by_state[state] = by_state.get(state, 0) + 1
        by_local[local] = by_local.get(local, 0) + 1
        by_cid_status[cid_status] = by_cid_status.get(cid_status, 0) + 1
        by_vendor_status[vendor_status] = by_vendor_status.get(vendor_status, 0) + 1
        if metrics.get("pubchem_cid_closed"):
            pubchem_cid_closed += 1
        if metrics.get("vendor_closed"):
            vendor_closed += 1
        if metrics.get("counts_as_organic_starting_material"):
            organic_starting_material_closures += 1
    return {
        "total_terminals": len(results),
        "pubchem_cid_closed": pubchem_cid_closed,
        "vendor_closed": vendor_closed,
        "organic_starting_material_closures": organic_starting_material_closures,
        "by_pubchem_cid_status": by_cid_status,
        "by_chemical_vendor_status": by_vendor_status,
        "by_buyability_state": by_state,
        "by_local_classification": by_local,
    }


def write_markdown(path: Path, payload: Dict[str, Any]) -> None:
    lines = [
        "# Terminal Buyability Audit",
        "",
        "This audit is conservative: PubChem CID evidence is not the same as live commercial availability, and metal-organic or inorganic reagent leaves are not treated as ordinary route terminals.",
        "",
        "## Summary",
        "",
    ]
    summary = payload.get("summary", {})
    lines.append(f"- total terminals: {summary.get('total_terminals', 0)}")
    lines.append(f"- PubChem CID closed: {summary.get('pubchem_cid_closed', 0)}")
    lines.append(f"- Vendor closed: {summary.get('vendor_closed', 0)}")
    for state, count in sorted((summary.get("by_buyability_state") or {}).items()):
        lines.append(f"- {state}: {count}")
    lines.extend(["", "## Terminals", ""])
    lines.append("| node | SMILES | local class | query | CID closed | Vendor closed | Evidence links | Allowlist evidence | Closure | decision |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for item in payload.get("results", []):
        local = item.get("local", {})
        pubchem = item.get("pubchem", {})
        metrics = item.get("pubchem_metrics", {})
        decision = item.get("buyability_decision", {})
        cids = ",".join(str(cid) for cid in pubchem.get("cids", [])[:3])
        query = pubchem.get("query_smiles", "") or "-"
        cid_cell = "yes" if metrics.get("pubchem_cid_closed") else "no"
        if metrics.get("cid_hit") and cids:
            cid_cell = f"yes ({cids})"
        vendor_cell = "yes" if metrics.get("vendor_closed") else "no"
        allowlist = item.get("allowlist", {})
        evidence = allowlist.get("evidence", {}) if isinstance(allowlist, dict) else {}
        link_parts = []
        cid_urls = metrics.get("cid_urls", []) if isinstance(metrics.get("cid_urls"), list) else []
        for cid, url in zip(pubchem.get("cids", [])[:3], cid_urls[:3]):
            if cid and url:
                link_parts.append(f"[CID {cid}]({url})")
        vendor_url = metrics.get("vendor_url", "")
        if vendor_url and (metrics.get("vendor_hit") or metrics.get("vendor_closed")):
            link_parts.append(f"[vendors]({vendor_url})")
        if allowlist.get("hit"):
            records = []
            if isinstance(evidence, dict):
                records.extend(evidence.get("forms", []) if isinstance(evidence.get("forms"), list) else [])
                records.extend(
                    evidence.get("source_reagents", [])
                    if isinstance(evidence.get("source_reagents"), list)
                    else []
                )
            for record in records:
                if not isinstance(record, dict):
                    continue
                cid = record.get("cid")
                cid_url = record.get("cid_url", "")
                if cid and cid_url:
                    link_parts.append(f"[{record.get('name', 'CID')} CID {cid}]({cid_url})")
                vendor = record.get("vendor_evidence", {})
                vendor_link = vendor.get("vendor_url", "") if isinstance(vendor, dict) else ""
                if record.get("vendor_hit") and vendor_link:
                    link_parts.append(f"[{record.get('name', 'vendor')}: vendors]({vendor_link})")
        links_cell = "<br/>".join(dict.fromkeys(link_parts[:6])) if link_parts else "-"
        allowlist_cell = evidence.get("basis", "") if allowlist.get("hit") else "-"
        closure_cell = metrics.get("closure_basis", "unknown")
        lines.append(
            "| {node} | `{smiles}` | {local} | `{query}` | {cid} | {vendor} | {links} | {allowlist} | {closure} | {decision} |".format(
                node=item.get("node_id", ""),
                smiles=item.get("smiles", ""),
                local=local.get("local_classification", ""),
                query=query,
                cid=cid_cell,
                vendor=vendor_cell,
                links=links_cell,
                allowlist=allowlist_cell,
                closure=closure_cell,
                decision=decision.get("state", ""),
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def default_output_paths(input_path: Path, terminal_path: Path) -> Tuple[Path, Path]:
    if input_path.is_dir():
        out_dir = input_path / "export" if (input_path / "export").exists() else input_path
    else:
        out_dir = terminal_path.parent
    return out_dir / "terminal_buyability_audit.json", out_dir / "terminal_buyability_audit.md"


def analysis_audit_output_paths(analysis_dir: Path, run_name: str) -> Tuple[Path, Path]:
    safe_name = "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in run_name)
    out_dir = analysis_dir / "run_audits"
    stem = f"{safe_name}_terminal_buyability_audit"
    return out_dir / f"{stem}.json", out_dir / f"{stem}.md"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Rachel export dir, terminals.json, or tree.json")
    parser.add_argument("--output-json", default="", help="Output JSON path")
    parser.add_argument("--output-md", default="", help="Output Markdown path")
    parser.add_argument("--cache-dir", default="", help="PubChem cache directory")
    parser.add_argument("--offline", action="store_true", help="Run local chemistry classification only")
    parser.add_argument("--no-vendors", action="store_true", help="Skip PUG View Chemical Vendors lookup")
    parser.add_argument("--query-reagents", action="store_true", help="Also query PubChem for reagent-classified leaves")
    parser.add_argument("--limit", type=int, default=0, help="Limit number of terminals for smoke tests")
    parser.add_argument("--timeout", type=float, default=20.0, help="PubChem request timeout seconds")
    parser.add_argument("--pause", type=float, default=0.2, help="Pause between PubChem requests")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    input_path = Path(args.input)
    records, terminal_path = load_terminal_records(input_path)
    if args.limit and args.limit > 0:
        records = records[: args.limit]

    default_json, default_md = default_output_paths(input_path, terminal_path)
    output_json = Path(args.output_json) if args.output_json else default_json
    output_md = Path(args.output_md) if args.output_md else default_md
    cache_dir = Path(args.cache_dir) if args.cache_dir else output_json.parent / ".pubchem_cache"

    client = PubChemClient(
        cache_dir=cache_dir,
        timeout=args.timeout,
        pause_seconds=args.pause,
        offline=args.offline,
    )
    results = [
        audit_record(
            record,
            client,
            include_vendors=not args.no_vendors and not args.offline,
            query_reagents=args.query_reagents,
        )
        for record in records
    ]
    payload = {
        "schema": "rachel-v2-terminal-buyability-audit-002",
        "input": str(terminal_path),
        "offline": bool(args.offline),
        "pubchem": {
            "pug_rest": f"{PUBCHEM_BASE}/pug",
            "pug_view": f"{PUBCHEM_BASE}/pug_view",
            "chemical_vendors_heading": "Chemical Vendors",
            "note": "CID/vendor evidence is not live inventory, price, or shipping confirmation.",
        },
        "summary": summarize(results),
        "results": results,
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    write_markdown(output_md, payload)
    print(json.dumps({
        "ok": True,
        "input": str(terminal_path),
        "output_json": str(output_json),
        "output_md": str(output_md),
        "summary": payload["summary"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
