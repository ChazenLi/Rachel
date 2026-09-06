"""Verify the extracted Rachel-beta-1.2-20260831_120000 release."""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_checksums() -> int:
    payload = json.loads((ROOT / "SHA256SUMS.json").read_text(encoding="utf-8"))
    for item in payload["files"]:
        path = ROOT / item["path"]
        if not path.is_file():
            raise SystemExit(f"checksum file missing: {item['path']}")
        if path.stat().st_size != item["bytes"]:
            raise SystemExit(f"size mismatch: {item['path']}")
        if sha256_file(path) != item["sha256"]:
            raise SystemExit(f"sha256 mismatch: {item['path']}")
    return len(payload["files"])


def main() -> None:
    checksum_count = verify_checksums()
    release_id = json.loads(
        (ROOT / "SHA256SUMS.json").read_text(encoding="utf-8")
    )["release"]

    from Rachel.main import RetroCmd
    from Rachel.main.prompt_mount import build_prompt_mount
    from Rachel.main.retro_cmd import PUBLIC_COMMANDS
    from Rachel.main.public_protocol import project_public_payload
    from Rachel.main.report_projection import build_route_report_view
    from Rachel.main.report_visualization import annotate_visualization
    from Rachel.chem_tools.atom_mapping_audit import audit_atom_mapping
    from Rachel.chem_tools.fg_delta_audit import audit_fg_delta
    from Rachel.chem_tools.graph_delta_audit import audit_graph_delta
    from Rachel.chem_tools.reaction_family_validate import validate_reaction_family
    from Rachel.chem_tools.ring_topology_audit import audit_ring_topology
    from Rachel.chem_tools.topology_intent import action_declares_topology_change
    from Rachel.chem_tools.validation_contract import build_validation_contract
    from Rachel.chem_tools.validation_policy import build_validation_gate
    from Rachel.tools.atom_bond_map import build_atom_bond_map
    from Rachel.tools.logged_runner import LoggedRetroCmd
    from Rachel.knowledge import get_base_profile
    from Rachel.knowledge.pack import KnowledgePackManager
    from Rachel.tools.audit_terminal_buyability_batch import build_arg_parser
    from Rachel.tools.pubchem_terminal_audit import load_terminal_allowlist
    from Rachel.tools.post_route_audit import POST_ROUTE_AUDIT_SCHEMA, route_digest

    validation_imports = [
        audit_atom_mapping,
        audit_fg_delta,
        audit_graph_delta,
        validate_reaction_family,
        audit_ring_topology,
        action_declares_topology_change,
        build_validation_contract,
        build_validation_gate,
        project_public_payload,
        build_route_report_view,
        annotate_visualization,
        KnowledgePackManager,
        route_digest,
        LoggedRetroCmd,
    ]
    if not all(callable(item) for item in validation_imports):
        raise SystemExit("validation dependency closure did not load")

    base_profile = get_base_profile()
    if len(base_profile.digest) != 64:
        raise SystemExit("base knowledge profile digest is invalid")
    if not base_profile.get("prompt.experience_cards").get("cards"):
        raise SystemExit("base knowledge profile did not load experience cards")
    prompt_mount = build_prompt_mount(
        "context_compact", knowledge_profile=base_profile
    )
    if prompt_mount.get("knowledge_profile_hash") != base_profile.digest:
        raise SystemExit("prompt mount did not retain the knowledge profile")
    if not prompt_mount.get("knowledge_refs"):
        raise SystemExit("prompt mount did not expose knowledge provenance")
    learning_commands = {
        "record_outcome",
        "learning_review",
        "propose_knowledge",
    }
    if not learning_commands.issubset(set(PUBLIC_COMMANDS)):
        raise SystemExit("controlled-learning commands are not public")
    if "audit" not in PUBLIC_COMMANDS:
        raise SystemExit("post-finalize audit command is not public")

    run = ROOT / "release_smoke"
    if run.exists():
        shutil.rmtree(run)
    run.mkdir()
    cmd = LoggedRetroCmd(str(run / "session.json"))

    init_result = cmd.execute("init", {
        "target": "CC(=O)NC",
        "name": "release_smoke_acetamide",
        "terminal_cs_threshold": 0.0,
        "max_depth": 4,
        "max_steps": 10,
    })
    if not init_result.get("ok"):
        raise SystemExit(f"init failed: {init_result}")

    next_result = cmd.execute("next", {})
    if "error" in next_result:
        raise SystemExit(f"next failed: {next_result}")

    structure_result = cmd.execute("context", {"detail": "structure"})
    structure_current = structure_result.get("current") or {}
    molecule_structure = structure_current.get("molecule_structure") or {}
    if len(molecule_structure.get("atoms", [])) <= 0:
        raise SystemExit(f"structure context missing atoms: {structure_result}")
    if len(molecule_structure.get("bonds", [])) <= 0:
        raise SystemExit(f"structure context missing bonds: {structure_result}")
    if "molecule_structure" in (next_result.get("current") or {}):
        raise SystemExit("compact context unexpectedly contains molecule_structure")

    atom_bond_map = build_atom_bond_map("CCO")
    if atom_bond_map.get("atom_count") != 3 or atom_bond_map.get("bond_count") != 2:
        raise SystemExit(f"atom_bond_map compatibility failed: {atom_bond_map}")

    plan_result = cmd.execute("route_plan", {
        "route_thesis": "Use a direct amide disconnection smoke strategy.",
        "route_mode": "late_fgi",
        "key_disconnections": ["amide C-N"],
        "mode_evidence": ["small amide smoke target"],
        "strategic_risks": ["smoke test only"],
        "revision_triggers": ["unexpected missing action-space"],
        "revision_reason": "release smoke",
    })
    if not plan_result.get("ok"):
        raise SystemExit(f"route_plan failed: {plan_result}")

    sites = cmd.execute("reaction_sites", {})
    if "error" in sites or int(sites.get("site_count", 0) or 0) <= 0:
        raise SystemExit(f"reaction_sites failed: {sites}")

    allowlist_count = len(load_terminal_allowlist())
    if allowlist_count <= 0:
        raise SystemExit("terminal allowlist did not load")
    args = build_arg_parser().parse_args(["--dataset", "n1", "--limit", "1"])
    if args.dataset != "n1" or args.limit != 1:
        raise SystemExit("terminal-audit parser verification failed")

    accept_result = cmd.execute("accept", {
        "reason": "Small stable material used only for release audit smoke.",
    })
    if not accept_result.get("ok"):
        raise SystemExit(f"release-smoke terminal acceptance failed: {accept_result}")
    queue_end = cmd.execute("next", {})
    if queue_end.get("action") != "queue_empty":
        raise SystemExit(f"release-smoke queue did not close: {queue_end}")
    finalized = cmd.execute("finalize", {"summary": "Offline post-route audit smoke."})
    if finalized.get("status", {}).get("status") != "complete":
        raise SystemExit(f"release-smoke finalize failed: {finalized}")

    digest_before_audit = route_digest(cmd.session.orch.tree)
    audit_result = cmd.execute("audit", {"offline": True})
    audit_payload = audit_result.get("post_route_audit") or {}
    if not audit_result.get("ok") or audit_payload.get("schema") != POST_ROUTE_AUDIT_SCHEMA:
        raise SystemExit(f"offline post-route audit failed: {audit_result}")
    audit_status = audit_result.get("post_route_audit_status") or {}
    if audit_status.get("state") != "current" or audit_status.get("audit_status") != "incomplete":
        raise SystemExit(f"offline post-route audit status is invalid: {audit_status}")
    if route_digest(cmd.session.orch.tree) != digest_before_audit:
        raise SystemExit("post-route audit changed route chemistry")

    export_dir = run / "export"
    export_result = cmd.execute("export", {
        "name": "release_smoke_acetamide",
        "output_dir": str(export_dir),
    })
    if export_result.get("post_route_audit_status", {}).get("state") != "current":
        raise SystemExit(f"current audit was not retained by export: {export_result}")
    audit_sidecars = [
        export_dir / "post_route_audit.json",
        export_dir / "post_route_audit.md",
    ]
    if not all(path.is_file() for path in audit_sidecars):
        raise SystemExit("current post-route audit sidecars were not exported")
    hidden_report_diagnostics = (
        "decision audit",
        "sandbox evidence",
        "prompt events",
        "experience cards",
        "custom provenance",
    )
    for report_name in (
        "report.txt",
        "SYNTHESIS_REPORT.md",
        "SYNTHESIS_REPORT.html",
    ):
        report_path = export_dir / report_name
        if not report_path.is_file():
            raise SystemExit(f"human report missing: {report_name}")
        rendered = report_path.read_text(encoding="utf-8").lower()
        if any(marker in rendered for marker in hidden_report_diagnostics):
            raise SystemExit(f"human report leaked diagnostics: {report_name}")
        if "post_route_audit.md" not in rendered:
            raise SystemExit(f"human report omitted audit sidecar link: {report_name}")
        if report_name.endswith(".html") and ("data:image" in rendered or "images/" not in rendered):
            raise SystemExit("HTML report did not reuse external image assets")

    reopened = cmd.execute("review_terminal", {
        "smiles": "CC(=O)NC",
        "reason": "Reopen only to verify route-digest staleness.",
    })
    if not reopened.get("ok"):
        raise SystemExit(f"release-smoke review_terminal failed: {reopened}")
    if cmd.execute("status", {}).get("post_route_audit", {}).get("state") != "stale":
        raise SystemExit("route mutation did not mark the post-route audit stale")
    stale_export_dir = run / "stale_export"
    stale_export = cmd.execute("export", {
        "name": "release_smoke_stale",
        "output_dir": str(stale_export_dir),
    })
    if stale_export.get("post_route_audit_status", {}).get("state") != "stale":
        raise SystemExit(f"stale audit status was not disclosed: {stale_export}")
    if any((stale_export_dir / name).exists() for name in ("post_route_audit.json", "post_route_audit.md")):
        raise SystemExit("stale post-route audit sidecars were exported")

    summary = {
        "ok": True,
        "release": release_id,
        "python": sys.executable,
        "checksum_files": checksum_count,
        "site_count": sites.get("site_count", 0),
        "structure_atoms": len(molecule_structure.get("atoms", [])),
        "structure_bonds": len(molecule_structure.get("bonds", [])),
        "atom_bond_map_atoms": atom_bond_map.get("atom_count", 0),
        "terminal_allowlist_entries": allowlist_count,
        "validation_modules": len(validation_imports),
        "knowledge_profile_hash": base_profile.digest,
        "controlled_learning_commands": sorted(learning_commands),
        "post_route_audit_command": "audit",
        "post_route_audit_schema": POST_ROUTE_AUDIT_SCHEMA,
        "post_route_audit_stale_export_blocked": True,
        "human_report_projection_verified": True,
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
