#!/usr/bin/env python
"""Print atom and bond index annotations for a SMILES string.

Examples:
    python E:/Python/skills/Rachel/tools/atom_bond_map.py --smiles "CCO"
    echo CCO | python E:/Python/skills/Rachel/tools/atom_bond_map.py
    python E:/Python/skills/Rachel/tools/atom_bond_map.py --smiles "CCO" --json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def build_atom_bond_map(smiles: str) -> Dict[str, Any]:
    """Return a structured atom/bond annotation for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")

    canonical_smiles = Chem.MolToSmiles(mol)
    atoms: List[Dict[str, Any]] = []
    bonds: List[Dict[str, Any]] = []

    for atom in mol.GetAtoms():
        atoms.append(
            {
                "idx": atom.GetIdx(),
                "symbol": atom.GetSymbol(),
                "is_aromatic": atom.GetIsAromatic(),
                "formal_charge": atom.GetFormalCharge(),
                "total_hs": atom.GetTotalNumHs(),
                "degree": atom.GetDegree(),
                "hybridization": str(atom.GetHybridization()),
            }
        )

    for bond in mol.GetBonds():
        bonds.append(
            {
                "idx": bond.GetIdx(),
                "begin_atom_idx": bond.GetBeginAtomIdx(),
                "end_atom_idx": bond.GetEndAtomIdx(),
                "bond_type": str(bond.GetBondType()),
                "is_aromatic": bond.GetIsAromatic(),
                "is_in_ring": bond.IsInRing(),
            }
        )

    return {
        "input_smiles": smiles,
        "canonical_smiles": canonical_smiles,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "atom_count": mol.GetNumAtoms(),
        "bond_count": mol.GetNumBonds(),
        "atoms": atoms,
        "bonds": bonds,
    }


def format_atom_bond_map(data: Dict[str, Any]) -> str:
    """Format the atom/bond annotation as plain text."""
    lines = [
        f"Input SMILES: {data['input_smiles']}",
        f"Canonical SMILES: {data['canonical_smiles']}",
        f"Formula: {data['formula']}",
        f"Atom count: {data['atom_count']}",
        f"Bond count: {data['bond_count']}",
        "---ATOMS---",
    ]

    for atom in data["atoms"]:
        aromatic = "arom" if atom["is_aromatic"] else ""
        lines.append(
            " ".join(
                [
                    str(atom["idx"]),
                    atom["symbol"],
                    aromatic,
                    "H",
                    str(atom["total_hs"]),
                    "charge",
                    str(atom["formal_charge"]),
                    "deg",
                    str(atom["degree"]),
                    atom["hybridization"],
                ]
            ).strip()
        )

    lines.append("---BONDS---")
    for bond in data["bonds"]:
        ring = "ring" if bond["is_in_ring"] else ""
        aromatic = "arom" if bond["is_aromatic"] else ""
        lines.append(
            " ".join(
                [
                    str(bond["idx"]),
                    str(bond["begin_atom_idx"]),
                    str(bond["end_atom_idx"]),
                    bond["bond_type"],
                    ring,
                    aromatic,
                ]
            ).strip()
        )

    return "\n".join(lines)


def _read_smiles(args: argparse.Namespace) -> str:
    if args.smiles:
        return args.smiles.strip()
    if not sys.stdin.isatty():
        return sys.stdin.read().strip()
    raise ValueError("provide --smiles or pipe a SMILES string on stdin")


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Print atom and bond index annotations for a SMILES string."
    )
    parser.add_argument(
        "--smiles",
        help="SMILES string to annotate. If omitted, the script reads from stdin.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit structured JSON instead of plain text.",
    )
    return parser.parse_args(argv)


def main(argv: List[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        smiles = _read_smiles(args)
        data = build_atom_bond_map(smiles)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if args.json:
        print(json.dumps(data, indent=2, ensure_ascii=False))
    else:
        print(format_atom_bond_map(data))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
