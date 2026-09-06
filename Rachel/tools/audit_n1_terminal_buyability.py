"""Compatibility wrapper for the dataset-parameterized terminal audit batch."""

from __future__ import annotations

import sys
from typing import Optional, List

from Rachel.tools.audit_terminal_buyability_batch import main as batch_main


def main(argv: Optional[List[str]] = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if "--dataset" not in args:
        args = ["--dataset", "n1", *args]
    return batch_main(args)


if __name__ == "__main__":
    raise SystemExit(main())
