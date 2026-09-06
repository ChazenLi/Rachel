"""Resolved Rachel knowledge profile with content verification and provenance."""

from __future__ import annotations

import copy
import hashlib
import json
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, Mapping, Optional


SCHEMA = "rachel.knowledge-pack.v1"
RESOURCE_SCHEMA = "rachel.knowledge-resource.v1"
_BASE_MANIFEST = Path(__file__).resolve().parent / "base" / "manifest.json"
_PACK_REF_RE = re.compile(r"^(?P<id>[A-Za-z0-9][A-Za-z0-9._-]*)@(?P<version>[A-Za-z0-9][A-Za-z0-9._-]*)$")


class KnowledgeProfileError(ValueError):
    """Raised when a knowledge pack is invalid or no longer reproducible."""


@dataclass(frozen=True)
class _PackIdentity:
    pack_id: str
    version: str
    kind: str


class KnowledgeProfile:
    """Read-only resolved knowledge with a small caller-facing interface."""

    def __init__(
        self,
        *,
        resources: Mapping[str, Any],
        sources: Mapping[str, Mapping[str, Dict[str, str]]],
        digest: str,
        pins: Optional[list[Mapping[str, str]]] = None,
    ) -> None:
        self._resources = copy.deepcopy(dict(resources))
        self._sources = copy.deepcopy(dict(sources))
        self.digest = digest
        self._pins = copy.deepcopy(list(pins or []))

    def get(self, resource: str) -> Any:
        """Return a defensive copy of one resolved resource."""
        if resource not in self._resources:
            raise KnowledgeProfileError(f"unknown knowledge resource: {resource}")
        return copy.deepcopy(self._resources[resource])

    def source(self, resource: str, entry_id: str) -> Dict[str, str]:
        """Return path-free provenance for one resolved entry."""
        try:
            return dict(self._sources[resource][entry_id])
        except KeyError as exc:
            raise KnowledgeProfileError(
                f"unknown knowledge entry: {resource}/{entry_id}"
            ) from exc

    def snapshot(self) -> Dict[str, Any]:
        """Return the path-free identity required to reproduce this profile."""
        return {
            "schema": "rachel.knowledge-profile.v1",
            "digest": self.digest,
            "packs": copy.deepcopy(self._pins),
        }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _entry_ids(data: Any, spec: Mapping[str, Any]) -> list[str]:
    merge = str(spec.get("merge", ""))
    if merge == "mapping":
        if not isinstance(data, dict):
            raise KnowledgeProfileError("mapping resource must contain a JSON object")
        return [str(key) for key in data if not str(key).startswith("__comment")]
    if merge == "catalog":
        path = str(spec.get("entries_path", ""))
        entries = data.get(path) if isinstance(data, dict) else None
        if not isinstance(entries, list):
            raise KnowledgeProfileError(f"catalog resource requires list field: {path}")
        id_field = str(spec.get("id_field", "id"))
        ids = [str(item.get(id_field, "")) for item in entries if isinstance(item, dict)]
        if not all(ids):
            raise KnowledgeProfileError("catalog resource contains an entry without an id")
        return ids
    if merge == "list_catalog":
        if not isinstance(data, list):
            raise KnowledgeProfileError("list catalog resource must contain a JSON list")
        id_field = str(spec.get("id_field", "id"))
        ids = [str(item.get(id_field, "")) for item in data if isinstance(item, dict)]
        if not all(ids):
            raise KnowledgeProfileError("list catalog contains an entry without an id")
        return ids
    raise KnowledgeProfileError(f"unsupported resource merge mode: {merge}")


def _load_base_profile(manifest_path: Path) -> KnowledgeProfile:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA:
        raise KnowledgeProfileError(f"unsupported knowledge manifest: {manifest.get('schema')}")
    identity = _PackIdentity(
        pack_id=str(manifest.get("id", "")),
        version=str(manifest.get("version", "")),
        kind=str(manifest.get("kind", "")),
    )
    if identity != _PackIdentity("rachel.base", "1.0.0", "base"):
        raise KnowledgeProfileError("built-in base manifest identity is invalid")

    resources: Dict[str, Any] = {}
    sources: Dict[str, Dict[str, Dict[str, str]]] = {}
    digest_resources: Dict[str, str] = {}
    for resource_name, spec in manifest.get("resources", {}).items():
        path = (manifest_path.parent / str(spec.get("path", ""))).resolve()
        expected_hash = str(spec.get("sha256", "")).lower()
        actual_hash = _sha256(path)
        if not expected_hash or actual_hash != expected_hash:
            raise KnowledgeProfileError(
                f"knowledge resource hash mismatch: {resource_name}"
            )
        data = json.loads(path.read_text(encoding="utf-8"))
        resources[str(resource_name)] = data
        provenance = {
            "pack_id": identity.pack_id,
            "pack_version": identity.version,
            "pack_kind": identity.kind,
            "resource": str(resource_name),
        }
        sources[str(resource_name)] = {
            entry_id: {**provenance, "entry_id": entry_id}
            for entry_id in _entry_ids(data, spec)
        }
        digest_resources[str(resource_name)] = actual_hash

    digest_payload = {
        "schema": SCHEMA,
        "packs": [{"id": identity.pack_id, "version": identity.version, "kind": identity.kind}],
        "resources": digest_resources,
    }
    digest = hashlib.sha256(
        json.dumps(
            digest_payload,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    return KnowledgeProfile(
        resources=resources,
        sources=sources,
        digest=digest,
        pins=[
            {
                "id": identity.pack_id,
                "version": identity.version,
                "kind": identity.kind,
                "manifest_sha256": _sha256(manifest_path),
            }
        ],
    )


def _parse_pack_ref(pack_ref: str) -> tuple[str, str]:
    match = _PACK_REF_RE.fullmatch(str(pack_ref).strip())
    if not match:
        raise KnowledgeProfileError(
            f"invalid knowledge pack reference (expected id@version): {pack_ref}"
        )
    return match.group("id"), match.group("version")


def _find_pack_manifest(pack_id: str, version: str, roots: list[Path]) -> Path:
    candidates: list[Path] = []
    for root in roots:
        for relative in (
            Path(pack_id) / version / "manifest.json",
            Path(f"{pack_id}@{version}") / "manifest.json",
        ):
            candidate = (root / relative).resolve()
            if candidate.is_file() and candidate not in candidates:
                candidates.append(candidate)
    if not candidates:
        raise KnowledgeProfileError(f"knowledge pack not found: {pack_id}@{version}")
    if len(candidates) > 1:
        raise KnowledgeProfileError(f"ambiguous knowledge pack: {pack_id}@{version}")
    return candidates[0]


def _load_external_resource(pack_dir: Path, spec: Mapping[str, Any], resource: str) -> Any:
    path = (pack_dir / str(spec.get("path", ""))).resolve()
    try:
        path.relative_to(pack_dir)
    except ValueError as exc:
        raise KnowledgeProfileError(
            f"knowledge resource escapes its pack directory: {resource}"
        ) from exc
    if path.suffix.lower() != ".json" or not path.is_file():
        raise KnowledgeProfileError(f"knowledge resource must be a JSON file: {resource}")
    expected_hash = str(spec.get("sha256", "")).lower()
    if not expected_hash or _sha256(path) != expected_hash:
        raise KnowledgeProfileError(f"knowledge resource hash mismatch: {resource}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema") != RESOURCE_SCHEMA:
        raise KnowledgeProfileError(f"invalid knowledge overlay schema: {resource}")
    if not isinstance(payload.get("entries"), list):
        raise KnowledgeProfileError(f"knowledge overlay requires entries: {resource}")
    return payload


def _merge_overlay(
    *,
    resource: str,
    base_spec: Mapping[str, Any],
    resolved: Any,
    resolved_sources: Dict[str, Dict[str, str]],
    overlay: Mapping[str, Any],
    identity: _PackIdentity,
) -> None:
    merge = str(base_spec.get("merge", ""))
    if merge == "catalog":
        entries_path = str(base_spec.get("entries_path", ""))
        target_entries = resolved.get(entries_path) if isinstance(resolved, dict) else None
        if not isinstance(target_entries, list):
            raise KnowledgeProfileError(f"invalid resolved catalog: {resource}")
        id_field = str(base_spec.get("id_field", "id"))
        existing = {
            str(item.get(id_field, ""))
            for item in target_entries
            if isinstance(item, dict)
        }
        positions = {
            str(item.get(id_field, "")): index
            for index, item in enumerate(target_entries)
            if isinstance(item, dict)
        }
    elif merge == "list_catalog":
        if not isinstance(resolved, list):
            raise KnowledgeProfileError(f"invalid resolved list catalog: {resource}")
        target_entries = resolved
        id_field = str(base_spec.get("id_field", "id"))
        existing = {
            str(item.get(id_field, ""))
            for item in target_entries
            if isinstance(item, dict)
        }
        positions = {
            str(item.get(id_field, "")): index
            for index, item in enumerate(target_entries)
            if isinstance(item, dict)
        }
    elif merge == "mapping":
        if not isinstance(resolved, dict):
            raise KnowledgeProfileError(f"invalid resolved mapping: {resource}")
        existing = set(resolved)
        positions = {}
    else:
        raise KnowledgeProfileError(f"unsupported resource merge mode: {merge}")

    namespace_prefixes = (f"{identity.pack_id}.", f"{identity.pack_id}:")
    for raw_entry in overlay.get("entries", []):
        if not isinstance(raw_entry, dict):
            raise KnowledgeProfileError(f"knowledge overlay entry must be an object: {resource}")
        metadata = raw_entry.get("_knowledge")
        if not isinstance(metadata, dict):
            raise KnowledgeProfileError(
                f"knowledge overlay entry requires explicit operation: {resource}"
            )
        operation = str(metadata.get("operation", ""))
        entry_id = str(raw_entry.get("id", "")).strip()
        clean_entry = {key: copy.deepcopy(value) for key, value in raw_entry.items() if key != "_knowledge"}

        if operation == "add":
            if not entry_id.startswith(namespace_prefixes):
                raise KnowledgeProfileError(
                    f"new knowledge entry must use pack namespace {identity.pack_id}: {entry_id}"
                )
            if entry_id in existing:
                raise KnowledgeProfileError(f"duplicate knowledge entry: {resource}/{entry_id}")
            if merge in {"catalog", "list_catalog"}:
                if id_field != "id":
                    clean_entry[id_field] = clean_entry.pop("id")
                positions[entry_id] = len(target_entries)
                target_entries.append(clean_entry)
            else:
                if "value" not in clean_entry:
                    raise KnowledgeProfileError(
                        f"mapping overlay entry requires value: {resource}/{entry_id}"
                    )
                resolved[entry_id] = clean_entry["value"]
            existing.add(entry_id)
        elif operation == "replace":
            target = str(metadata.get("target", "")).strip()
            reason = str(metadata.get("reason", "")).strip()
            supersedes = str(metadata.get("supersedes", "")).strip()
            if not target or not reason or supersedes != target:
                raise KnowledgeProfileError(
                    f"replace requires target, matching supersedes, and reason: {resource}"
                )
            if target not in existing or entry_id != target:
                raise KnowledgeProfileError(f"replace target not found: {resource}/{target}")
            if bool(base_spec.get("locked")):
                raise KnowledgeProfileError(f"locked knowledge resource cannot be replaced: {resource}")
            if merge in {"catalog", "list_catalog"}:
                current = target_entries[positions[target]]
                if isinstance(current, dict) and current.get("locked") is True:
                    raise KnowledgeProfileError(
                        f"locked knowledge entry cannot be replaced: {resource}/{target}"
                    )
                if id_field != "id":
                    clean_entry[id_field] = clean_entry.pop("id")
                target_entries[positions[target]] = clean_entry
            else:
                current = resolved[target]
                if isinstance(current, dict) and current.get("locked") is True:
                    raise KnowledgeProfileError(
                        f"locked knowledge entry cannot be replaced: {resource}/{target}"
                    )
                if "value" not in clean_entry:
                    raise KnowledgeProfileError(
                        f"mapping replacement requires value: {resource}/{target}"
                    )
                resolved[target] = clean_entry["value"]
            entry_id = target
        elif operation == "disable":
            target = str(metadata.get("target", "")).strip()
            reason = str(metadata.get("reason", "")).strip()
            supersedes = str(metadata.get("supersedes", "")).strip()
            if not target or not reason or supersedes != target:
                raise KnowledgeProfileError(
                    f"disable requires target, matching supersedes, and reason: {resource}"
                )
            if target not in existing or entry_id != target:
                raise KnowledgeProfileError(f"disable target not found: {resource}/{target}")
            if bool(base_spec.get("locked")):
                raise KnowledgeProfileError(f"locked knowledge resource cannot be disabled: {resource}")
            if merge in {"catalog", "list_catalog"}:
                current = target_entries[positions[target]]
                if isinstance(current, dict) and current.get("locked") is True:
                    raise KnowledgeProfileError(
                        f"locked knowledge entry cannot be disabled: {resource}/{target}"
                    )
                target_entries.pop(positions[target])
                positions = {
                    str(item.get(id_field, "")): index
                    for index, item in enumerate(target_entries)
                    if isinstance(item, dict)
                }
            else:
                current = resolved[target]
                if isinstance(current, dict) and current.get("locked") is True:
                    raise KnowledgeProfileError(
                        f"locked knowledge entry cannot be disabled: {resource}/{target}"
                    )
                resolved.pop(target)
            existing.remove(target)
            resolved_sources.pop(target, None)
            continue
        else:
            raise KnowledgeProfileError(
                f"unsupported knowledge operation: {resource}/{operation or '<missing>'}"
            )

        resolved_sources[entry_id] = {
            "pack_id": identity.pack_id,
            "pack_version": identity.version,
            "pack_kind": identity.kind,
            "resource": resource,
            "entry_id": entry_id,
        }


def resolve_knowledge_profile(
    pack_refs: Optional[list[str]] = None,
    *,
    knowledge_roots: Optional[list[Path | str]] = None,
) -> KnowledgeProfile:
    """Resolve built-in base plus ordered versioned team/project packs."""
    refs = list(pack_refs or [])
    if not refs:
        return get_base_profile()
    roots = [Path(root).resolve() for root in knowledge_roots or []]
    if not roots:
        raise KnowledgeProfileError("knowledge_roots required for external packs")

    base = get_base_profile()
    resources = copy.deepcopy(base._resources)
    sources = copy.deepcopy(base._sources)
    base_manifest = json.loads(_BASE_MANIFEST.read_text(encoding="utf-8"))
    base_specs = base_manifest.get("resources", {})
    loaded = {"rachel.base@1.0.0"}
    pack_digests: list[Dict[str, Any]] = []
    project_seen = False

    for pack_ref in refs:
        pack_id, version = _parse_pack_ref(pack_ref)
        manifest_path = _find_pack_manifest(pack_id, version, roots)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("schema") != SCHEMA:
            raise KnowledgeProfileError(f"unsupported knowledge manifest: {pack_ref}")
        kind = str(manifest.get("kind", ""))
        if manifest.get("id") != pack_id or manifest.get("version") != version:
            raise KnowledgeProfileError(f"knowledge pack identity mismatch: {pack_ref}")
        if kind not in {"team", "project"}:
            raise KnowledgeProfileError(f"external knowledge pack kind must be team or project: {pack_ref}")
        if project_seen and kind == "team":
            raise KnowledgeProfileError("team knowledge packs must precede project packs")
        project_seen = project_seen or kind == "project"
        dependencies = manifest.get("dependencies", [])
        if not isinstance(dependencies, list) or any(dep not in loaded for dep in dependencies):
            raise KnowledgeProfileError(f"unresolved knowledge pack dependency: {pack_ref}")
        identity = _PackIdentity(pack_id, version, kind)
        resource_digests: Dict[str, str] = {}
        for resource, spec in manifest.get("resources", {}).items():
            if resource not in base_specs:
                raise KnowledgeProfileError(f"unknown knowledge resource: {resource}")
            overlay = _load_external_resource(manifest_path.parent, spec, resource)
            _merge_overlay(
                resource=resource,
                base_spec=base_specs[resource],
                resolved=resources[resource],
                resolved_sources=sources[resource],
                overlay=overlay,
                identity=identity,
            )
            resource_digests[str(resource)] = str(spec.get("sha256", "")).lower()
        manifest_hash = _sha256(manifest_path)
        pack_digests.append(
            {
                "id": pack_id,
                "version": version,
                "kind": kind,
                "manifest_sha256": manifest_hash,
                "resources": resource_digests,
            }
        )
        loaded.add(f"{pack_id}@{version}")

    digest_payload = {
        "schema": SCHEMA,
        "base_digest": base.digest,
        "packs": pack_digests,
    }
    digest = hashlib.sha256(
        json.dumps(digest_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return KnowledgeProfile(
        resources=resources,
        sources=sources,
        digest=digest,
        pins=base._pins + [
            {
                "id": str(item["id"]),
                "version": str(item["version"]),
                "kind": str(item["kind"]),
                "manifest_sha256": str(item["manifest_sha256"]),
            }
            for item in pack_digests
        ],
    )


def resolve_pinned_knowledge_profile(
    snapshot: Mapping[str, Any],
    *,
    knowledge_roots: Optional[list[Path | str]] = None,
) -> KnowledgeProfile:
    """Resolve and verify a previously pinned profile without hot updates."""
    if snapshot.get("schema") != "rachel.knowledge-profile.v1":
        raise KnowledgeProfileError("unsupported pinned knowledge profile")
    pins = snapshot.get("packs")
    if not isinstance(pins, list) or not pins:
        raise KnowledgeProfileError("pinned knowledge profile has no packs")
    base_pin = pins[0] if isinstance(pins[0], dict) else {}
    if (
        base_pin.get("id") != "rachel.base"
        or base_pin.get("version") != "1.0.0"
        or base_pin.get("kind") != "base"
    ):
        raise KnowledgeProfileError("pinned knowledge profile is missing base")
    refs = [
        f"{pin.get('id')}@{pin.get('version')}"
        for pin in pins[1:]
        if isinstance(pin, dict)
    ]
    profile = resolve_knowledge_profile(refs, knowledge_roots=knowledge_roots)
    actual = profile.snapshot()
    if actual.get("digest") != snapshot.get("digest") or actual.get("packs") != pins:
        raise KnowledgeProfileError("pinned knowledge profile hash mismatch")
    return profile


@lru_cache(maxsize=1)
def get_base_profile() -> KnowledgeProfile:
    """Return the verified built-in base profile."""
    return _load_base_profile(_BASE_MANIFEST)
