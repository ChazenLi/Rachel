"""Batch terminal buyability audit for Rachel-beta de novo runs.

This writes all per-run audit JSON/Markdown files under the analysis directory,
not under individual walkthrough run exports.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from Rachel.tools.pubchem_terminal_audit import (
    PUBCHEM_BASE,
    PubChemClient,
    analysis_audit_output_paths,
    audit_record,
    load_terminal_records,
    summarize,
    write_markdown,
)


def _bool_text(value: Any) -> str:
    return "true" if bool(value) else "false"


def _route_id_sort_key(value: str) -> tuple[int, str]:
    return (int(value), value) if value.isdigit() else (10**9, value)


def build_run_re(dataset: str) -> re.Pattern[str]:
    escaped = re.escape(dataset.lower())
    return re.compile(rf"^\d{{8}}(?:_\d{{6}})?_{escaped}_(\d+)(?:_de_novo)?$")


def find_run_dirs(walkthrough_runs: Path, dataset: str) -> List[Path]:
    run_re = build_run_re(dataset)
    return sorted(
        [
            item
            for item in walkthrough_runs.iterdir()
            if item.is_dir() and run_re.match(item.name.lower())
        ],
        key=lambda path: path.name,
    )


def expected_route_ids(token_csv: Path, dataset: str) -> List[str]:
    if not token_csv.exists():
        return []
    route_ids: List[str] = []
    with token_csv.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if str(row.get("dataset", "")).strip().lower() != dataset.lower():
                continue
            route_id = str(row.get("route_id") or row.get("id") or "").strip()
            if route_id:
                route_ids.append(route_id)
    return sorted(set(route_ids), key=_route_id_sort_key)


def _write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_batch(
    *,
    dataset: str,
    walkthrough_runs: Path,
    analysis_dir: Path,
    token_csv: Path,
    cache_dir: Path,
    timeout: float,
    pause: float,
    limit: int = 0,
) -> Dict[str, Any]:
    client = PubChemClient(
        cache_dir=cache_dir,
        timeout=timeout,
        pause_seconds=pause,
        offline=False,
    )
    dataset = dataset.lower()
    run_re = build_run_re(dataset)
    run_dirs = find_run_dirs(walkthrough_runs, dataset)
    if limit > 0:
        run_dirs = run_dirs[:limit]
    run_rows: List[Dict[str, Any]] = []
    terminal_rows: List[Dict[str, Any]] = []

    for index, run_dir in enumerate(run_dirs, 1):
        match = run_re.match(run_dir.name.lower())
        if match is None:
            continue
        route_id = match.group(1)
        records, terminal_path = load_terminal_records(run_dir)
        results = [
            audit_record(record, client, include_vendors=True, query_reagents=False)
            for record in records
        ]
        payload = {
            "schema": "rachel-v2-terminal-buyability-audit-002",
            "input": str(terminal_path),
            "run": run_dir.name,
            "route_id": route_id,
            "offline": False,
            "pubchem": {
                "pug_rest": f"{PUBCHEM_BASE}/pug",
                "pug_view": f"{PUBCHEM_BASE}/pug_view",
                "chemical_vendors_heading": "Chemical Vendors",
                "note": "CID/vendor evidence is not live inventory, price, or shipping confirmation.",
            },
            "summary": summarize(results),
            "results": results,
        }
        output_json, output_md = analysis_audit_output_paths(analysis_dir, run_dir.name)
        output_json.parent.mkdir(parents=True, exist_ok=True)
        output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        write_markdown(output_md, payload)

        summary = payload["summary"]
        total = int(summary.get("total_terminals", 0))
        cid_closed = int(summary.get("pubchem_cid_closed", 0))
        vendor_closed = int(summary.get("vendor_closed", 0))
        run_rows.append({
            "run": run_dir.name,
            "route_id": route_id,
            "total_terminals": total,
            "pubchem_cid_closed_terminals": cid_closed,
            "vendor_closed_terminals": vendor_closed,
            "pubchem_cid_closed_run": _bool_text(cid_closed == total),
            "vendor_closed_run": _bool_text(vendor_closed == total),
            "audit_json": str(output_json),
            "audit_md": str(output_md),
            "status": "audited",
        })

        for item in results:
            metrics = item.get("pubchem_metrics", {})
            pubchem = item.get("pubchem", {})
            decision = item.get("buyability_decision", {})
            terminal_rows.append({
                "run": run_dir.name,
                "route_id": route_id,
                "node_id": item.get("node_id", ""),
                "smiles": item.get("smiles", ""),
                "local_classification": item.get("local", {}).get("local_classification", ""),
                "query_smiles": pubchem.get("query_smiles", ""),
                "pubchem_cid_closed": _bool_text(metrics.get("pubchem_cid_closed")),
                "vendor_closed": _bool_text(metrics.get("vendor_closed")),
                "closure_basis": metrics.get("closure_basis", ""),
                "decision": decision.get("state", ""),
                "closure_evidence_json": json.dumps(
                    metrics.get("closure_evidence", {}),
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
            })
        print(
            f"[{index}/{len(run_dirs)}] {run_dir.name}: "
            f"CID {cid_closed}/{total}, vendor {vendor_closed}/{total}"
        )

    expected_ids = expected_route_ids(token_csv, dataset)
    audited_ids = sorted({row["route_id"] for row in run_rows}, key=_route_id_sort_key)
    missing_ids = [route_id for route_id in expected_ids if route_id not in set(audited_ids)]

    route_rows = list(run_rows)
    for route_id in missing_ids:
        route_rows.append({
            "run": "",
            "route_id": route_id,
            "total_terminals": "",
            "pubchem_cid_closed_terminals": "",
            "vendor_closed_terminals": "",
            "pubchem_cid_closed_run": "",
            "vendor_closed_run": "",
            "audit_json": "",
            "audit_md": "",
            "status": "missing_export",
        })
    route_rows.sort(key=lambda row: (row["status"] == "missing_export", _route_id_sort_key(str(row["route_id"]))))

    total_terminals = len(terminal_rows)
    cid_closed_terminals = sum(row["pubchem_cid_closed"] == "true" for row in terminal_rows)
    vendor_closed_terminals = sum(row["vendor_closed"] == "true" for row in terminal_rows)
    cid_closed_runs = sum(row["pubchem_cid_closed_run"] == "true" for row in run_rows)
    vendor_closed_runs = sum(row["vendor_closed_run"] == "true" for row in run_rows)
    not_cid_rows = [row for row in terminal_rows if row["pubchem_cid_closed"] != "true"]
    not_vendor_rows = [row for row in terminal_rows if row["vendor_closed"] != "true"]
    not_vendor_runs = [row for row in run_rows if row["vendor_closed_run"] != "true"]

    run_fields = [
        "run", "route_id", "total_terminals", "pubchem_cid_closed_terminals",
        "vendor_closed_terminals", "pubchem_cid_closed_run", "vendor_closed_run",
        "audit_json", "audit_md", "status",
    ]
    terminal_fields = [
        "run", "route_id", "node_id", "smiles", "local_classification",
        "query_smiles", "pubchem_cid_closed", "vendor_closed", "closure_basis",
        "decision", "closure_evidence_json",
    ]
    _write_csv(analysis_dir / f"{dataset}_terminal_buyability_runs.csv", run_rows, run_fields)
    _write_csv(analysis_dir / f"{dataset}_terminal_buyability_route_ids.csv", route_rows, run_fields)
    _write_csv(analysis_dir / f"{dataset}_terminal_buyability_terminals.csv", terminal_rows, terminal_fields)

    summary = {
        "schema": "rachel-v2-terminal-buyability-batch-summary-001",
        "dataset": dataset,
        "export_runs_audited": len(run_rows),
        "unique_export_route_ids": len(set(audited_ids)),
        "token_csv_route_ids": len(expected_ids),
        "missing_token_route_ids": missing_ids,
        "terminals_audited": total_terminals,
        "pubchem_cid_closed_terminals": cid_closed_terminals,
        "vendor_closed_terminals": vendor_closed_terminals,
        "pubchem_cid_closed_runs": cid_closed_runs,
        "vendor_closed_runs": vendor_closed_runs,
        "not_pubchem_cid_closed_terminals": not_cid_rows,
        "not_vendor_closed_terminals": not_vendor_rows,
        "not_vendor_closed_runs": not_vendor_runs,
    }
    (analysis_dir / f"{dataset}_terminal_buyability_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    lines = [
        f"# {dataset.upper()} Terminal Buyability Summary",
        "",
        "- closure metrics: `pubchem_cid_closed` and `vendor_closed`",
        "- allowlist entries are evidence mappings inside `closure_evidence`, not a third closure metric",
        "",
        f"- export runs audited: {len(run_rows)}",
        f"- unique export route IDs: {len(set(audited_ids))}",
        f"- token.csv {dataset.upper()} route IDs: {len(expected_ids)}",
        f"- token IDs missing export: {', '.join(missing_ids) if missing_ids else 'none'}",
        f"- terminals audited: {total_terminals}",
        f"- PubChem CID closed terminals: {cid_closed_terminals}/{total_terminals}",
        f"- Vendor closed terminals: {vendor_closed_terminals}/{total_terminals}",
        f"- PubChem CID closed runs: {cid_closed_runs}/{len(run_rows)}",
        f"- Vendor closed runs: {vendor_closed_runs}/{len(run_rows)}",
        "",
        "## Runs Not Vendor Closed",
        "",
        "| run | terminals | CID closed | Vendor closed |",
        "|---|---:|---:|---:|",
    ]
    for row in not_vendor_runs:
        lines.append(
            f"| {row['run']} | {row['total_terminals']} | "
            f"{row['pubchem_cid_closed_terminals']} | {row['vendor_closed_terminals']} |"
        )
    lines.extend([
        "",
        "## Terminal Rows Not Vendor Closed",
        "",
        "| run | node | SMILES | CID closed | Vendor closed | decision |",
        "|---|---|---|---|---|---|",
    ])
    for row in not_vendor_rows:
        lines.append(
            f"| {row['run']} | {row['node_id']} | `{row['smiles']}` | "
            f"{row['pubchem_cid_closed']} | {row['vendor_closed']} | {row['decision']} |"
        )
    if not_cid_rows:
        lines.extend([
            "",
            "## Terminal Rows Not PubChem CID Closed",
            "",
            "| run | node | SMILES | decision |",
            "|---|---|---|---|",
        ])
        for row in not_cid_rows:
            lines.append(f"| {row['run']} | {row['node_id']} | `{row['smiles']}` | {row['decision']} |")
    (analysis_dir / f"{dataset}_terminal_buyability_summary.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    parser.add_argument("--dataset", default="n1", help="Dataset label, for example n1 or n5")
    parser.add_argument("--walkthrough-runs", default=str(root / "walkthrough_runs"))
    parser.add_argument("--analysis-dir", default="")
    parser.add_argument("--token-csv", default="")
    parser.add_argument("--cache-dir", default="")
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument("--pause", type=float, default=0.05)
    parser.add_argument("--limit", type=int, default=0)
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    root = Path(__file__).resolve().parents[2]
    dataset = str(args.dataset).strip().lower()
    if not dataset:
        raise ValueError("--dataset must not be empty")
    walkthrough_runs = Path(args.walkthrough_runs)
    token_csv = Path(args.token_csv) if args.token_csv else walkthrough_runs / "token.csv"
    analysis_dir = Path(args.analysis_dir) if args.analysis_dir else root / "analysis" / f"{dataset}_terminal_buyability"
    cache_dir = Path(args.cache_dir) if args.cache_dir else root / "analysis" / "pubchem_terminal_audit_cache"
    started = time.time()
    summary = run_batch(
        dataset=dataset,
        walkthrough_runs=walkthrough_runs,
        analysis_dir=analysis_dir,
        token_csv=token_csv,
        cache_dir=cache_dir,
        timeout=args.timeout,
        pause=args.pause,
        limit=args.limit,
    )
    print(json.dumps({
            "ok": True,
            "dataset": dataset,
            "elapsed_sec": round(time.time() - started, 1),
            "summary": {
                "export_runs_audited": summary["export_runs_audited"],
                "missing_token_route_ids": summary["missing_token_route_ids"],
            "terminals_audited": summary["terminals_audited"],
            "pubchem_cid_closed_terminals": summary["pubchem_cid_closed_terminals"],
            "vendor_closed_terminals": summary["vendor_closed_terminals"],
            "pubchem_cid_closed_runs": summary["pubchem_cid_closed_runs"],
            "vendor_closed_runs": summary["vendor_closed_runs"],
        },
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
