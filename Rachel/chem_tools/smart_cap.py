"""
Smart Capping — 智能断键推理
============================
当模板不命中或 LLM 需要辅助时，根据键两端的化学环境自动推断
合理的 capping 基团，生成前体 SMILES。

核心函数:
  - suggest_capping(smiles, bond) → List[CappingProposal]

设计原则:
  - 基于键两端原子类型、杂化、邻近官能团的规则匹配
  - 不依赖模板库，是独立的化学推理层
  - 输出带置信度的多种 capping 方案，供 LLM 选择或修改
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import RWMol

from ._rdkit_utils import parse_mol, canonical
from Rachel.knowledge.conditions import condition_matches


@dataclass
class CappingProposal:
    """一种 capping 方案。"""
    reaction_type: str          # e.g. "Suzuki coupling"
    fragments: List[str]        # 两端的前体 SMILES
    confidence: float           # 0-1
    description: str = ""       # 简短说明
    cap_atoms: Dict[int, str] = field(default_factory=dict)  # {atom_idx: cap_group}

    def to_dict(self) -> Dict[str, Any]:
        return {
            "reaction_type": self.reaction_type,
            "fragments": self.fragments,
            "confidence": round(self.confidence, 3),
            "description": self.description,
        }


# ─────────────────────────────────────────────────────────────────────────
# 键环境分析
# ─────────────────────────────────────────────────────────────────────────

def _atom_env(mol: Chem.Mol, idx: int) -> Dict[str, Any]:
    """分析单个原子的化学环境。"""
    atom = mol.GetAtomWithIdx(idx)
    symbol = atom.GetSymbol()
    hyb = str(atom.GetHybridization())
    in_ring = atom.IsInRing()
    aromatic = atom.GetIsAromatic()
    degree = atom.GetDegree()
    n_hs = atom.GetTotalNumHs()
    formal_charge = atom.GetFormalCharge()

    # 邻居原子信息
    neighbors = []
    for nb in atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(idx, nb.GetIdx())
        neighbors.append({
            "idx": nb.GetIdx(),
            "symbol": nb.GetSymbol(),
            "aromatic": nb.GetIsAromatic(),
            "bond_type": str(bond.GetBondType()) if bond else "",
        })

    # 检测特殊环境
    is_carbonyl_c = False
    is_carbamate_c = False
    is_amide_n = False
    is_ester_o = False
    is_carbamate_o = False
    is_carboxylic_acid_o = False
    adjacent_carbonyl = False

    if symbol == "C" and not aromatic:
        for nb in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(idx, nb.GetIdx())
            if nb.GetSymbol() == "O" and bond and "DOUBLE" in str(bond.GetBondType()):
                is_carbonyl_c = True
                break
        if is_carbonyl_c:
            has_single_bond_n = any(
                nb.GetSymbol() == "N"
                and mol.GetBondBetweenAtoms(idx, nb.GetIdx()).GetBondType()
                == Chem.BondType.SINGLE
                for nb in atom.GetNeighbors()
            )
            has_single_bond_o = any(
                nb.GetSymbol() == "O"
                and mol.GetBondBetweenAtoms(idx, nb.GetIdx()).GetBondType()
                == Chem.BondType.SINGLE
                for nb in atom.GetNeighbors()
            )
            is_carbamate_c = has_single_bond_n and has_single_bond_o

    if symbol == "N":
        for nb in atom.GetNeighbors():
            nb_atom = mol.GetAtomWithIdx(nb.GetIdx())
            if nb_atom.GetSymbol() == "C":
                for nb2 in nb_atom.GetNeighbors():
                    bond2 = mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx())
                    if nb2.GetSymbol() == "O" and bond2 and "DOUBLE" in str(bond2.GetBondType()):
                        is_amide_n = True
                        break

    if symbol == "O" and not aromatic:
        for nb in atom.GetNeighbors():
            nb_atom = mol.GetAtomWithIdx(nb.GetIdx())
            if nb_atom.GetSymbol() == "C":
                for nb2 in nb_atom.GetNeighbors():
                    bond2 = mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx())
                    if nb2.GetSymbol() == "O" and bond2 and "DOUBLE" in str(bond2.GetBondType()):
                        is_ester_o = True
                        is_carbamate_o = any(
                            other.GetSymbol() == "N"
                            and mol.GetBondBetweenAtoms(
                                nb.GetIdx(), other.GetIdx()
                            ).GetBondType() == Chem.BondType.SINGLE
                            for other in nb_atom.GetNeighbors()
                        )
                        is_carboxylic_acid_o = (
                            degree == 1
                            and (n_hs > 0 or formal_charge < 0)
                            and not is_carbamate_o
                        )
                        if is_carboxylic_acid_o:
                            is_ester_o = False
                        break

    # 相邻是否有 C=O
    for nb in atom.GetNeighbors():
        nb_atom = mol.GetAtomWithIdx(nb.GetIdx())
        if nb_atom.GetSymbol() == "C":
            for nb2 in nb_atom.GetNeighbors():
                bond2 = mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx())
                if nb2.GetSymbol() == "O" and bond2 and "DOUBLE" in str(bond2.GetBondType()):
                    adjacent_carbonyl = True

    return {
        "idx": idx,
        "symbol": symbol,
        "hybridization": hyb,
        "in_ring": in_ring,
        "aromatic": aromatic,
        "degree": degree,
        "n_hs": n_hs,
        "formal_charge": formal_charge,
        "neighbors": neighbors,
        "is_carbonyl_c": is_carbonyl_c,
        "is_carbamate_c": is_carbamate_c,
        "is_amide_n": is_amide_n,
        "is_ester_o": is_ester_o,
        "is_carbamate_o": is_carbamate_o,
        "is_carboxylic_acid_o": is_carboxylic_acid_o,
        "adjacent_carbonyl": adjacent_carbonyl,
        "is_sp2": "SP2" in hyb,
        "is_sp3": "SP3" in hyb,
        "has_halide_neighbor": any(
            item["symbol"] in {"Br", "Cl", "I"} for item in neighbors
        ),
    }


def _bond_env(mol: Chem.Mol, i: int, j: int) -> Dict[str, Any]:
    """分析键的化学环境。"""
    bond = mol.GetBondBetweenAtoms(i, j)
    if bond is None:
        return {"error": f"no bond between {i} and {j}"}

    bt = str(bond.GetBondType())
    in_ring = bond.IsInRing()
    conjugated = bond.GetIsConjugated()

    env_i = _atom_env(mol, i)
    env_j = _atom_env(mol, j)

    return {
        "bond_type": bt,
        "in_ring": in_ring,
        "conjugated": conjugated,
        "atom_i": env_i,
        "atom_j": env_j,
    }


# ─────────────────────────────────────────────────────────────────────────
# Capping 规则引擎
# ─────────────────────────────────────────────────────────────────────────


def _sanitize_capped(rw) -> Optional[str]:
    """对 capping 后的 RWMol 做多级 sanitize，返回 canonical SMILES 或 None。

    Level 1: 标准 SanitizeMol
    Level 2: 跳过 SANITIZE_PROPERTIES
    Level 3: 清除残留芳香标记后重试（处理断开芳香环后的片段）
    """
    from rdkit.Chem import RWMol as _RW

    # Level 1
    try:
        Chem.SanitizeMol(rw)
        return Chem.MolToSmiles(rw, canonical=True)
    except Exception:
        pass

    # Level 2
    try:
        Chem.SanitizeMol(
            rw,
            Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
        )
        return Chem.MolToSmiles(rw, canonical=True)
    except Exception:
        pass

    # Level 3: 清除非环原子的芳香标记
    try:
        rw2 = _RW(rw.GetMol()) if not isinstance(rw, _RW) else _RW(rw)
        ri = rw2.GetRingInfo()
        ring_atoms = set()
        for ring in ri.AtomRings():
            ring_atoms.update(ring)
        for atom in rw2.GetAtoms():
            if atom.GetIdx() not in ring_atoms and atom.GetIsAromatic():
                atom.SetIsAromatic(False)
        for bond in rw2.GetBonds():
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a1 not in ring_atoms or a2 not in ring_atoms:
                if bond.GetBondType() == Chem.BondType.AROMATIC:
                    bond.SetBondType(Chem.BondType.SINGLE)
        Chem.SanitizeMol(rw2)
        return Chem.MolToSmiles(rw2, canonical=True)
    except Exception:
        pass

    # Level 4: 全部清除芳香标记，让 RDKit 重新推断
    try:
        rw3 = _RW(rw.GetMol()) if not isinstance(rw, _RW) else _RW(rw)
        for atom in rw3.GetAtoms():
            atom.SetIsAromatic(False)
        for bond in rw3.GetBonds():
            if bond.GetBondType() == Chem.BondType.AROMATIC:
                bond.SetBondType(Chem.BondType.SINGLE)
        Chem.SanitizeMol(rw3)
        return Chem.MolToSmiles(rw3, canonical=True)
    except Exception:
        return None


def _replace_dummy_with_cap(frag_input, cap_smiles: str,
                            target_isotope: int = 0) -> Optional[str]:
    """将片段中的 dummy atom 替换为 capping 基团（纯图操作）。

    Parameters
    ----------
    frag_input : str | Chem.Mol
        含 dummy atom 的片段。可以是 SMILES 字符串或 RDKit Mol 对象。
        对于芳香环内键断开的片段，必须传 Mol 对象（SMILES round-trip 会失败）。
    cap_smiles : str
        要加的基团。支持:
        - "[H]" / "H"  — 加氢
        - "=O"          — 加双键氧（如还原胺化逆向 N-CH₂R → NH + R'CHO）
        - 单原子: "Br", "Cl", "I", "F", "O", "N", "S" 等
        - 多原子 SMILES: "B(O)O", "[Mg]Br", "[Zn]Cl", "OC(=O)" 等
    target_isotope : int
        如果 > 0，只替换该 isotope 的 dummy（用于环内键场景，
        同一片段有两个 dummy 需要分别替换）。0 = 替换所有 dummy。
    """
    try:
        # 支持 str 或 Mol 输入
        if isinstance(frag_input, Chem.Mol):
            frag_mol = Chem.RWMol(frag_input)
            input_is_mol = True
        else:
            frag_mol = Chem.MolFromSmiles(frag_input)
            input_is_mol = False
        if frag_mol is None:
            return None

        # 找到目标 dummy atom（原子序数 0）
        dummy_indices = []
        for atom in frag_mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                if target_isotope == 0 or atom.GetIsotope() == target_isotope:
                    dummy_indices.append(atom.GetIdx())

        if not dummy_indices:
            # 没有匹配的 dummy，原样返回
            if input_is_mol:
                return _sanitize_capped(Chem.RWMol(frag_mol))
            return frag_input

        from rdkit.Chem import RWMol

        # ── 特殊处理: [H] — 加氢 ──
        if cap_smiles in ("[H]", "H"):
            rw = RWMol(frag_mol)
            for d_idx in sorted(dummy_indices, reverse=True):
                dummy_atom = rw.GetAtomWithIdx(d_idx)
                for nb in dummy_atom.GetNeighbors():
                    nb.SetNumExplicitHs(nb.GetNumExplicitHs() + 1)
                rw.RemoveAtom(d_idx)
            return _sanitize_capped(rw)

        # ── 特殊处理: =O — 加双键氧 ──
        if cap_smiles == "=O":
            rw = RWMol(frag_mol)
            for d_idx in sorted(dummy_indices, reverse=True):
                dummy_atom = rw.GetAtomWithIdx(d_idx)
                nb_list = list(dummy_atom.GetNeighbors())
                if not nb_list:
                    return None
                nb_idx = nb_list[0].GetIdx()
                rw.RemoveAtom(d_idx)
                new_nb_idx = nb_idx if nb_idx < d_idx else nb_idx - 1
                o_idx = rw.AddAtom(Chem.Atom(8))
                rw.AddBond(new_nb_idx, o_idx, Chem.BondType.DOUBLE)
            return _sanitize_capped(rw)

        # ── 通用路径: 图操作替换 dummy → cap 分子 ──
        cap_mol = Chem.MolFromSmiles(cap_smiles)
        if cap_mol is None:
            return None

        rw = RWMol(frag_mol)

        for d_idx in sorted(dummy_indices, reverse=True):
            dummy_atom = rw.GetAtomWithIdx(d_idx)
            nb_list = list(dummy_atom.GetNeighbors())
            if not nb_list:
                return None
            # dummy 只连一个邻居（断键端原子）
            anchor_idx = nb_list[0].GetIdx()
            # dummy-anchor 之间的键类型（通常是 SINGLE）
            d_bond = rw.GetBondBetweenAtoms(d_idx, anchor_idx)
            bond_type = d_bond.GetBondType() if d_bond else Chem.BondType.SINGLE

            # 删除 dummy
            rw.RemoveAtom(d_idx)
            # 删除后 anchor_idx 可能需要调整
            new_anchor = anchor_idx if anchor_idx < d_idx else anchor_idx - 1

            if cap_mol.GetNumAtoms() == 1:
                # 单原子 cap: 直接加一个原子
                cap_atom = cap_mol.GetAtomWithIdx(0)
                new_atom = Chem.Atom(cap_atom.GetAtomicNum())
                new_atom.SetFormalCharge(cap_atom.GetFormalCharge())
                new_idx = rw.AddAtom(new_atom)
                rw.AddBond(new_anchor, new_idx, bond_type)
            else:
                # 多原子 cap: 合并 cap 分子，连接到 anchor
                # cap 的连接点 = 第一个原子（atom index 0）
                offset = rw.GetNumAtoms()
                combined = Chem.CombineMols(rw.GetMol(), cap_mol)
                rw = RWMol(combined)
                cap_connect = offset  # cap 分子的第一个原子在合并后的索引
                rw.AddBond(new_anchor, cap_connect, bond_type)

        return _sanitize_capped(rw)

    except Exception:
        return None




# ─────────────────────────────────────────────────────────────────────────
# Capping 规则表
# ─────────────────────────────────────────────────────────────────────────

# 每条规则: (条件函数, reaction_type, cap_for_i, cap_for_j, confidence, description)
# 条件函数接收 (env_i, env_j, bond_env) 返回 bool
# cap_for_i: 给 i 端加的基团 SMILES
# cap_for_j: 给 j 端加的基团 SMILES

def _build_profile_capping_rules(
    env: Dict[str, Any], knowledge_profile=None
) -> List[Dict[str, Any]]:
    """Evaluate resolved declarative rules against one bond environment."""
    if knowledge_profile is None:
        from Rachel.knowledge import get_base_profile

        knowledge_profile = get_base_profile()
    resource = knowledge_profile.get("chem.smart_cap_rules")
    entries = [
        entry
        for entry in resource.get("entries", []) or []
        if isinstance(entry, dict) and entry.get("id")
    ]
    entries.sort(
        key=lambda entry: (
            -int(entry.get("priority", 0) or 0),
            str(entry.get("id", "")),
        )
    )

    matched = [
        entry
        for entry in entries
        if entry.get("fallback_rank") is None
        and condition_matches(entry.get("condition"), env)
    ]
    if not matched:
        fallbacks = [
            entry
            for entry in entries
            if entry.get("fallback_rank") is not None
            and condition_matches(entry.get("condition"), env)
        ]
        if fallbacks:
            first_rank = min(int(entry["fallback_rank"]) for entry in fallbacks)
            matched = [
                entry
                for entry in fallbacks
                if int(entry["fallback_rank"]) == first_rank
            ]

    proposals: List[Dict[str, Any]] = []
    for entry in matched:
        source = knowledge_profile.source(
            "chem.smart_cap_rules", str(entry["id"])
        )
        candidates = entry.get("proposals") or [entry]
        for candidate in candidates:
            if not isinstance(candidate, dict):
                continue
            proposals.append(
                {
                    "reaction_type": str(candidate.get("reaction_type", "")),
                    "cap_i": str(candidate.get("cap_i", "")),
                    "cap_j": str(candidate.get("cap_j", "")),
                    "confidence": float(candidate.get("confidence", 0.0) or 0.0),
                    "description": str(candidate.get("description", "")),
                    "knowledge_ref": source,
                }
            )
    return proposals


def suggest_capping(
    smiles: str,
    bond: Tuple[int, int],
    max_proposals: int = 5,
    knowledge_profile=None,
) -> Dict[str, Any]:
    """为指定键位推断合理的 capping 方案。

    根据键两端原子的化学环境（元素、杂化、芳香性、邻近官能团），
    应用规则引擎推断断键后两端应该加什么基团。

    Parameters
    ----------
    smiles : str
        目标分子 SMILES。
    bond : Tuple[int, int]
        要断的键的原子索引对 (i, j)。
    max_proposals : int
        最多返回几个方案。

    Returns
    -------
    dict
        {
            "ok": True,
            "smiles": str,
            "bond": [i, j],
            "bond_env": {...},  # 键环境分析
            "proposals": [
                {
                    "reaction_type": str,
                    "fragments": [str, str],  # 两个前体 SMILES
                    "confidence": float,
                    "description": str,
                },
                ...
            ],
            "n_proposals": int,
        }
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": f"invalid SMILES: {smiles}"}

    i, j = int(bond[0]), int(bond[1])
    n_atoms = mol.GetNumAtoms()
    if i < 0 or j < 0 or i >= n_atoms or j >= n_atoms:
        return {"ok": False, "error": f"atom index out of range: ({i}, {j})"}

    bond_obj = mol.GetBondBetweenAtoms(i, j)
    if bond_obj is None:
        return {"ok": False, "error": f"no bond between atoms {i} and {j}"}

    # 分析键环境
    env = _bond_env(mol, i, j)

    # 生成 capping 规则
    rules = _build_profile_capping_rules(env, knowledge_profile)

    # 对每条规则，尝试生成前体 SMILES
    proposals: List[Dict[str, Any]] = []

    for rule in rules:
        cap_i = rule["cap_i"]
        cap_j = rule["cap_j"]

        # 断键并 cap
        frags = _generate_capped_fragments(smiles, i, j, cap_i, cap_j)
        if frags is None:
            continue

        # 验证所有前体 SMILES 合法性
        all_valid = all(parse_mol(f) is not None for f in frags)
        if not all_valid:
            continue

        proposals.append({
            "reaction_type": rule["reaction_type"],
            "fragments": [canonical(f) or f for f in frags],
            "confidence": rule["confidence"],
            "description": rule["description"],
            "ring_opening": len(frags) == 1,
            "knowledge_ref": rule["knowledge_ref"],
        })

    # 去重（按 fragments 集合）
    seen = set()
    unique_proposals = []
    for p in proposals:
        key = frozenset(p["fragments"])
        if key not in seen:
            seen.add(key)
            unique_proposals.append(p)

    # 按 confidence 排序
    unique_proposals.sort(key=lambda x: -x["confidence"])
    unique_proposals = unique_proposals[:max_proposals]

    return {
        "ok": True,
        "smiles": canonical(smiles) or smiles,
        "bond": [i, j],
        "bond_env": {
            "atom_i": {"symbol": env["atom_i"]["symbol"],
                       "aromatic": env["atom_i"]["aromatic"],
                       "hybridization": env["atom_i"]["hybridization"],
                       "is_carbonyl_c": env["atom_i"].get("is_carbonyl_c", False)},
            "atom_j": {"symbol": env["atom_j"]["symbol"],
                       "aromatic": env["atom_j"]["aromatic"],
                       "hybridization": env["atom_j"]["hybridization"],
                       "is_carbonyl_c": env["atom_j"].get("is_carbonyl_c", False)},
            "bond_type": env["bond_type"],
            "in_ring": env["in_ring"],
        },
        "proposals": unique_proposals,
        "n_proposals": len(unique_proposals),
    }


def custom_cap(
    smiles: str,
    bond: Tuple[int, int],
    cap_i: str,
    cap_j: str,
    reaction_type: str = "llm_custom",
) -> Dict[str, Any]:
    """LLM 自定义 capping：指定两端 cap 基团，系统执行断键+替换+验证。

    LLM 知道该加什么基团但不想手算完整前体 SMILES 时使用。
    比 try_precursors 轻量（只需指定 cap），比 suggest_capping 灵活（不受规则限制）。

    Parameters
    ----------
    smiles : str
        目标分子 SMILES。
    bond : Tuple[int, int]
        要断的键 (atom_i, atom_j)。
    cap_i : str
        给 atom_i 端加的基团 SMILES，如 "Br", "O", "B(O)O", "[H]", "=O", "Cl" 等。
    cap_j : str
        给 atom_j 端加的基团 SMILES。
    reaction_type : str
        反应类型标签（仅用于记录）。

    Returns
    -------
    dict
        {
            "ok": True,
            "fragments": [frag_i, frag_j],
            "cap_i": str, "cap_j": str,
            "reaction_type": str,
        }
    """
    mol = parse_mol(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid SMILES: %s" % smiles}

    i, j = int(bond[0]), int(bond[1])
    n_atoms = mol.GetNumAtoms()
    if i < 0 or j < 0 or i >= n_atoms or j >= n_atoms:
        return {"ok": False, "error": "atom index out of range: (%d, %d)" % (i, j)}

    bond_obj = mol.GetBondBetweenAtoms(i, j)
    if bond_obj is None:
        return {"ok": False, "error": "no bond between atoms %d and %d" % (i, j)}

    frags = _generate_capped_fragments(smiles, i, j, cap_i, cap_j)
    if frags is None:
        return {"ok": False, "error": "capping failed (cap_i=%s, cap_j=%s)" % (cap_i, cap_j)}

    # 验证所有前体合法性
    for idx, f in enumerate(frags):
        if parse_mol(f) is None:
            return {"ok": False, "error": "invalid fragment_%d: %s" % (idx, f)}

    can_frags = [canonical(f) or f for f in frags]
    is_ring_opening = len(frags) == 1

    result = {
        "ok": True,
        "smiles": canonical(smiles) or smiles,
        "bond": [i, j],
        "cap_i": cap_i,
        "cap_j": cap_j,
        "fragments": can_frags,
        "reaction_type": reaction_type,
        "ring_opening": is_ring_opening,
    }
    return result


def _generate_capped_fragments(
    smiles: str, i: int, j: int,
    cap_i: str, cap_j: str,
) -> Optional[Tuple[str, ...]]:
    """断键后对两端分别加 capping 基团。

    使用 FragmentOnBonds 断键（生成 dummy atoms），
    然后将 dummy 替换为 capping 基团。

    cap_i: 给 atom i 端加的基团
    cap_j: 给 atom j 端加的基团

    返回:
      - 非环键: (frag_i_smi, frag_j_smi) — 两个前体
      - 环内键: (single_smi,) — 一个前体（环打开后两端分别 cap）

    FragmentOnBonds dummy 标签逻辑:
      dummyLabels=(1, 2) 中 1=begin label, 2=end label
      begin/end 是 RDKit 内部的键方向，不一定等于 i/j
      关键: 含 [1*] dummy 的片段 = 不含 begin atom 的片段 = 含 end atom 的片段
            含 [2*] dummy 的片段 = 不含 end atom 的片段 = 含 begin atom 的片段
    """
    mol = parse_mol(smiles)
    if mol is None:
        return None

    bond = mol.GetBondBetweenAtoms(i, j)
    if bond is None:
        return None

    try:
        bond_idx = bond.GetIdx()
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        fragmented = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True,
                                          dummyLabels=[(1, 2)])
        if fragmented is None:
            return None

        frags = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=False)

        # ── 环内键：1 个片段，2 个 dummy ──
        if len(frags) == 1:
            frag = frags[0]
            dummies = {}
            for atom in frag.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    dummies[atom.GetIsotope()] = atom.GetIdx()
            if len(dummies) < 2:
                return None

            # isotope 2 dummy 邻接 begin atom，isotope 1 dummy 邻接 end atom
            # 确定哪个 isotope 对应 atom i / atom j
            if begin_idx == i:
                # begin=i → isotope 2 dummy 在 i 端, isotope 1 dummy 在 j 端
                iso_for_i = 2
                iso_for_j = 1
            else:
                iso_for_i = 1
                iso_for_j = 2

            # 直接传 Mol 对象，避免芳香环片段 SMILES round-trip 失败
            result = _replace_dummy_with_cap(frag, cap_i, target_isotope=iso_for_i)
            if result is None:
                return None
            result = _replace_dummy_with_cap(result, cap_j, target_isotope=iso_for_j)
            if result is None:
                return None

            return (result,)

        # ── 非环键：2 个片段 ──
        if len(frags) < 2:
            return None

        # 确定哪个片段包含 atom i，哪个包含 atom j
        # isotope==2 dummy → 这个片段包含 begin atom
        # isotope==1 dummy → 这个片段包含 end atom
        frag_with_begin = None
        frag_with_end = None

        for frag in frags:
            for atom in frag.GetAtoms():
                if atom.GetAtomicNum() == 0:  # dummy
                    isotope = atom.GetIsotope()
                    if isotope == 2:
                        frag_with_begin = frag
                    elif isotope == 1:
                        frag_with_end = frag

        if frag_with_begin is None or frag_with_end is None:
            # fallback
            frag_with_begin = frags[0]
            frag_with_end = frags[1]

        # 映射 begin/end → i/j
        if begin_idx == i:
            frag_with_i = frag_with_begin
            frag_with_j = frag_with_end
        else:
            frag_with_i = frag_with_end
            frag_with_j = frag_with_begin

        # 将 dummy 替换为 cap
        smi_i = Chem.MolToSmiles(frag_with_i)
        smi_j = Chem.MolToSmiles(frag_with_j)

        capped_i = _replace_dummy_with_cap(smi_i, cap_i)
        capped_j = _replace_dummy_with_cap(smi_j, cap_j)

        if capped_i is None or capped_j is None:
            return None

        return (capped_i, capped_j)

    except Exception:
        return None


