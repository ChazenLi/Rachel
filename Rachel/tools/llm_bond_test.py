#!/usr/bin/env python
"""LLM 交互断键测试平台 — 手动输入 SMILES 并测试断键操作的命令行工具。

用法:
    python Rachel/tools/llm_bond_test.py

交互流程:
    1. 输入产物 SMILES
    2. 查看所有可断键位及反应模板
    3. 选择键索引对 (如 "3 7") 或反应类型名称
    4. 查看断键结果
    5. 重复或输入 quit/exit 退出
"""

from __future__ import annotations

import sys
import os

# Ensure the project root is on sys.path so Rachel package is importable
_project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from Rachel.chem_tools.template_scan import find_disconnectable_bonds
from Rachel.chem_tools.bond_break import execute_disconnection
from Rachel.chem_tools._rdkit_utils import parse_mol, canonical


# ── Formatting helpers ─────────────────────────────────────────────────────

def _print_header(text: str) -> None:
    width = max(len(text) + 4, 50)
    print("\n" + "=" * width)
    print(f"  {text}")
    print("=" * width)


def _print_bond_table(bonds: list) -> None:
    """Pretty-print the list of disconnectable bonds."""
    print(f"\n  {'#':<4} {'Atoms':<10} {'Type':<8} {'Ring':<6} {'Score':<8} {'Templates'}")
    print("  " + "-" * 70)
    for idx, bond in enumerate(bonds):
        atoms = bond.get("atoms", [])
        btype = bond.get("bond_type", "?")
        in_ring = "Yes" if bond.get("in_ring") else "No"
        score = bond.get("heuristic_score", 0.0)
        templates = bond.get("templates", [])
        tpl_names = ", ".join(t.get("name", t.get("template_id", "?")) for t in templates)
        if not tpl_names:
            tpl_names = "(no template)"
        print(f"  {idx:<4} {str(atoms):<10} {btype:<8} {in_ring:<6} {score:<8.4f} {tpl_names}")


def _print_templates_for_bond(bond: dict) -> None:
    """Print detailed template info for a single bond."""
    templates = bond.get("templates", [])
    if not templates:
        print("  No templates associated with this bond.")
        return
    for i, tpl in enumerate(templates):
        print(f"\n  Template {i}:")
        print(f"    Name:       {tpl.get('name', '?')}")
        print(f"    Category:   {tpl.get('category', '?')}")
        print(f"    Confidence: {tpl.get('confidence', '?')}")
        if tpl.get("conditions"):
            print(f"    Conditions: {tpl['conditions']}")


def _print_break_result(result) -> None:
    """Print the BreakResult."""
    if result.success:
        print("\n  ✓ Disconnection succeeded!")
        print(f"  Precursors: {' . '.join(result.precursors)}")
    else:
        print("\n  ✗ Disconnection failed.")
        print(f"  Reason: {result.error}")
    if result.operation_log:
        print("\n  Operation log:")
        for entry in result.operation_log:
            print(f"    - {entry}")


# ── Bond selection parsing ─────────────────────────────────────────────────

def _parse_bond_selection(user_input: str, bonds: list):
    """Parse user input to select a bond.

    Accepts:
      - A pair of atom indices like "3 7" or "3,7"
      - A bond table index like "#0" or just "0" (if purely numeric)
      - A reaction type name (fuzzy match against template names)

    Returns (atoms_tuple, reaction_type) or None on failure.
    """
    text = user_input.strip()

    # Try "#N" table index
    if text.startswith("#"):
        try:
            idx = int(text[1:])
            if 0 <= idx < len(bonds):
                bond = bonds[idx]
                atoms = tuple(bond["atoms"])
                tpls = bond.get("templates", [])
                rtype = tpls[0].get("name", "") if tpls else ""
                return atoms, rtype
        except ValueError:
            pass

    # Try "i j" or "i,j" atom pair
    parts = text.replace(",", " ").split()
    if len(parts) == 2:
        try:
            i, j = int(parts[0]), int(parts[1])
            # Find matching bond entry
            for bond in bonds:
                a = bond["atoms"]
                if (a[0] == i and a[1] == j) or (a[0] == j and a[1] == i):
                    tpls = bond.get("templates", [])
                    rtype = tpls[0].get("name", "") if tpls else ""
                    return (i, j), rtype
            # Bond not in discovered list — still allow it
            return (i, j), ""
        except ValueError:
            pass

    # Try single number as table index
    if text.isdigit():
        idx = int(text)
        if 0 <= idx < len(bonds):
            bond = bonds[idx]
            atoms = tuple(bond["atoms"])
            tpls = bond.get("templates", [])
            rtype = tpls[0].get("name", "") if tpls else ""
            return atoms, rtype

    # Try reaction type name — find first bond with matching template
    key = text.lower()
    for bond in bonds:
        for tpl in bond.get("templates", []):
            name = tpl.get("name", "").lower()
            cat = tpl.get("category", "").lower()
            tid = tpl.get("template_id", "").lower()
            if key in name or key in cat or key in tid:
                return tuple(bond["atoms"]), tpl.get("name", "")
    return None


# ── Main interactive loop ──────────────────────────────────────────────────

def main() -> None:
    print()
    _print_header("LLM Bond Disconnection Test Platform")
    print("""
  Interactive tool for testing bond disconnection on molecules.

  Commands:
    - Enter a SMILES string to analyze
    - Type 'quit' or 'exit' to leave
    - Press Ctrl+C at any time to exit
""")

    while True:
        # ── Step 1: Get SMILES ──
        try:
            smiles_input = input("  SMILES> ").strip()
        except (KeyboardInterrupt, EOFError):
            print("\n\n  Goodbye!")
            break

        if not smiles_input:
            continue
        if smiles_input.lower() in ("quit", "exit", "q"):
            print("\n  Goodbye!")
            break

        # Validate SMILES
        mol = parse_mol(smiles_input)
        if mol is None:
            print(f"\n  ✗ Invalid SMILES: '{smiles_input}'")
            print("  Please enter a valid SMILES string.\n")
            continue

        can_smi = canonical(smiles_input) or smiles_input
        print(f"\n  Canonical SMILES: {can_smi}")
        print(f"  Heavy atoms: {mol.GetNumHeavyAtoms()}")

        # ── Step 2: Find disconnectable bonds ──
        result = find_disconnectable_bonds(smiles_input)
        if not result.get("ok"):
            print(f"\n  ✗ find_disconnectable_bonds failed: {result.get('error', 'unknown')}\n")
            continue

        bonds = result.get("bonds", [])
        unmatched = result.get("unmatched_bonds", [])

        if not bonds and not unmatched:
            print("\n  No disconnectable bonds found for this molecule.\n")
            continue

        if bonds:
            _print_header("Disconnectable Bonds (with templates)")
            _print_bond_table(bonds)

        if unmatched:
            print(f"\n  + {len(unmatched)} strategic bonds without template matches")

        # ── Step 3: Bond selection loop ──
        while True:
            print("\n  Select a bond to disconnect:")
            print("    - Enter atom pair (e.g. '3 7')")
            print("    - Enter table index (e.g. '0' or '#0')")
            print("    - Enter reaction type name (e.g. 'suzuki')")
            print("    - Enter 'back' for a new SMILES")
            print("    - Enter 'detail N' to see templates for bond #N")

            try:
                sel_input = input("  Bond> ").strip()
            except (KeyboardInterrupt, EOFError):
                print("\n\n  Goodbye!")
                return

            if not sel_input:
                continue
            if sel_input.lower() in ("back", "b"):
                print()
                break
            if sel_input.lower() in ("quit", "exit", "q"):
                print("\n  Goodbye!")
                return

            # Handle "detail N"
            if sel_input.lower().startswith("detail"):
                parts = sel_input.split()
                if len(parts) >= 2:
                    try:
                        didx = int(parts[1])
                        if 0 <= didx < len(bonds):
                            _print_templates_for_bond(bonds[didx])
                        else:
                            print(f"  Invalid index. Range: 0-{len(bonds)-1}")
                    except ValueError:
                        print("  Usage: detail <index>")
                else:
                    print("  Usage: detail <index>")
                continue

            # Parse selection
            parsed = _parse_bond_selection(sel_input, bonds)
            if parsed is None:
                print(f"  ✗ Could not parse bond selection: '{sel_input}'")
                print(f"  Valid indices: 0-{len(bonds)-1}")
                continue

            bond_pair, reaction_type = parsed
            print(f"\n  Executing disconnection: bond={bond_pair}, type='{reaction_type}'")

            # ── Step 4: Execute disconnection ──
            try:
                break_result = execute_disconnection(
                    smiles_input, bond_pair, reaction_type
                )
                _print_break_result(break_result)
            except Exception as exc:
                print(f"\n  ✗ Exception during disconnection: {exc}")

            print()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n  Interrupted. Goodbye!")
        sys.exit(0)
