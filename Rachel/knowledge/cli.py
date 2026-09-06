"""Manage Rachel knowledge packs as immutable, auditable JSON assets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from .pack import KnowledgePackError, KnowledgePackManager


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("validate", help="validate one published pack")
    validate.add_argument("--pack", type=Path, required=True)
    diff = commands.add_parser("diff", help="compare two published packs")
    diff.add_argument("--left", type=Path, required=True)
    diff.add_argument("--right", type=Path, required=True)
    conflicts = commands.add_parser(
        "conflicts", help="find blockers in an ordered pack composition"
    )
    conflicts.add_argument("--pack", type=Path, action="append", required=True)
    import_drafts = commands.add_parser(
        "import-drafts", help="import inactive session drafts into staging"
    )
    import_drafts.add_argument("--session", type=Path, required=True)
    import_drafts.add_argument("--workspace", type=Path, required=True)
    import_drafts.add_argument("--pack-id", required=True)
    for command in ("approve", "reject"):
        review = commands.add_parser(command, help=f"{command} one staging draft")
        review.add_argument("--workspace", type=Path, required=True)
        review.add_argument("--draft-id", required=True)
        review.add_argument("--reviewer", required=True)
        review.add_argument("--reason", required=True)
    publish = commands.add_parser("publish", help="publish a new immutable pack version")
    publish.add_argument("--workspace", type=Path, required=True)
    publish.add_argument("--pack-root", type=Path, required=True)
    publish.add_argument("--version", required=True)
    publish.add_argument("--kind", choices=("team", "project"), required=True)
    publish.add_argument("--dependency", action="append", default=[])
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    manager = KnowledgePackManager()
    try:
        if args.command == "validate":
            payload = manager.validate(args.pack)
        elif args.command == "diff":
            payload = manager.diff(args.left, args.right)
        elif args.command == "conflicts":
            payload = manager.conflicts(args.pack)
        elif args.command == "import-drafts":
            payload = manager.import_drafts(
                session=args.session,
                workspace=args.workspace,
                pack_id=args.pack_id,
            )
        elif args.command in {"approve", "reject"}:
            payload = manager.review_draft(
                workspace=args.workspace,
                draft_id=args.draft_id,
                decision="approved" if args.command == "approve" else "rejected",
                reviewer=args.reviewer,
                reason=args.reason,
            )
        elif args.command == "publish":
            payload = manager.publish(
                workspace=args.workspace,
                pack_root=args.pack_root,
                version=args.version,
                kind=args.kind,
                dependencies=args.dependency,
            )
        else:  # pragma: no cover - argparse prevents this path
            raise KnowledgePackError(f"unsupported command: {args.command}")
    except KnowledgePackError as exc:
        print(json.dumps({"ok": False, "error": str(exc)}, ensure_ascii=False))
        return 1
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0 if payload.get("ok") else 1


if __name__ == "__main__":
    raise SystemExit(main())
