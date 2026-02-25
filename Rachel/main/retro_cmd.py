#!/usr/bin/env python
"""
retro_cmd.py — LLM 多步逆合成命令接口
======================================
标准化的 JSON-in / JSON-out 命令层。
LLM 通过 Python 调用此脚本的函数，所有 SMILES 和上下文
只在 JSON 文件中流转，不经过终端。

命令列表:
  init          创建新会话
  next          取下一个分子（自动跳过 quick_pass）
  context       获取当前上下文（compact/full/status/tree）
  explore       展开键位详情
  explore_fgi   展开 FGI 详情
  try_bond      沙盒试断键
  try_fgi       沙盒试 FGI
  try_precursors 沙盒试 LLM 自提前体
  sandbox_list  查看沙盒历史
  sandbox_clear 清空沙盒
  select        选中沙盒方案
  commit        提交决策写入树
  accept        标记当前分子为 terminal
  skip          跳过当前分子
  tree          打印合成树
  status        查看编排状态
  finalize      完成编排
  report        生成正向合成报告

用法:
  python -m Rachel.main.retro_cmd <command> [session_file] [args_json]

  或在 Python 中:
    from Rachel.main.retro_cmd import RetroCmd
    cmd = RetroCmd("session.json")
    result = cmd.execute("init", {"target": "CC(=O)Nc1ccc(O)cc1", "name": "Paracetamol"})
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# 确保项目根目录在 path 中
ROOT = Path(__file__).resolve().parent.parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Rachel.main.retro_session import RetroSession
from Rachel.main.retro_report import generate_forward_report, get_terminal_list
from Rachel.main.retro_output import export_results


class RetroCmd:
    """LLM 逆合成命令接口。

    所有方法返回 dict，可直接 json.dumps。
    session_file 是唯一的状态载体。
    """

    def __init__(self, session_file: str = ".rachel_session.json"):
        self.session_file = session_file
        self._session: Optional[RetroSession] = None

    @property
    def session(self) -> RetroSession:
        if self._session is None:
            self._session = RetroSession.load(self.session_file)
        return self._session

    def execute(self, command: str, args: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """统一入口。command → 对应方法。"""
        args = args or {}
        dispatch = {
            "init": self.cmd_init,
            "next": self.cmd_next,
            "context": self.cmd_context,
            "explore": self.cmd_explore,
            "explore_fgi": self.cmd_explore_fgi,
            "try_bond": self.cmd_try_bond,
            "try_fgi": self.cmd_try_fgi,
            "try_precursors": self.cmd_try_precursors,
            "sandbox_list": self.cmd_sandbox_list,
            "sandbox_clear": self.cmd_sandbox_clear,
            "select": self.cmd_select,
            "commit": self.cmd_commit,
            "accept": self.cmd_accept,
            "skip": self.cmd_skip,
            "tree": self.cmd_tree,
            "status": self.cmd_status,
            "finalize": self.cmd_finalize,
            "report": self.cmd_report,
            "export": self.cmd_export,
            "smart_cap": self.cmd_smart_cap,
            "custom_cap": self.cmd_custom_cap,
        }
        fn = dispatch.get(command)
        if fn is None:
            return {"error": f"unknown command: {command}", "available": list(dispatch.keys())}
        try:
            return fn(args)
        except Exception as e:
            return {"error": str(e), "command": command}

    # ── 会话管理 ──

    def cmd_init(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """创建新会话。
        args: target (SMILES), name, max_depth, max_steps, terminal_cs_threshold
        """
        target = args.get("target", "")
        if not target:
            return {"error": "target SMILES required"}

        self._session = RetroSession.create(
            target,
            session_path=self.session_file,
            target_name=args.get("name", ""),
            max_depth=args.get("max_depth", 15),
            max_steps=args.get("max_steps", 50),
            terminal_cs_threshold=args.get("terminal_cs_threshold", 2.5),
        )
        return {
            "ok": True,
            "session_id": self._session.session_id,
            "session_file": str(self._session.path),
            "target": target,
            "name": args.get("name", ""),
            "hint": "调用 next 获取第一个分子的上下文",
        }

    # ── 编排流程 ──

    def cmd_next(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """取下一个分子。自动跳过所有 quick_pass terminal，
        直到遇到 standard 分子或队列为空。

        返回 compact 上下文 + 沙盒状态。
        """
        auto_terminals = []
        while True:
            ctx = self.session.prepare_next()
            if ctx is None:
                if auto_terminals:
                    return {
                        "action": "queue_empty",
                        "auto_terminals": auto_terminals,
                        "hint": "队列已空。调用 finalize 完成编排，或 tree 查看当前树。",
                    }
                return {
                    "action": "queue_empty",
                    "hint": "队列已空。调用 finalize 完成编排。",
                }

            if ctx.get("action") == "auto_terminal":
                auto_terminals.append({
                    "smiles": ctx["smiles"],
                    "cs_score": ctx["cs_score"],
                })
                continue

            # standard 分子
            result = ctx
            if auto_terminals:
                result["auto_terminals_skipped"] = auto_terminals
            return result

    def cmd_context(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """获取当前上下文。detail: compact/full/status/tree"""
        detail = args.get("detail", "compact")
        return self.session.get_context(detail)

    def cmd_explore(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """展开键位详情。args: bond_idx"""
        bond_idx = args.get("bond_idx", 0)
        return self.session.explore_bond(bond_idx)

    def cmd_explore_fgi(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """展开 FGI 详情。"""
        return self.session.explore_fgi()

    # ── 沙盒 ──

    def cmd_try_bond(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """沙盒试断键。args: bond_idx, alt_idx"""
        return self.session.sandbox_try(
            bond_idx=args.get("bond_idx", 0),
            alt_idx=args.get("alt_idx", 0),
        )

    def cmd_try_fgi(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """沙盒试 FGI。args: fgi_idx"""
        return self.session.sandbox_try(fgi_idx=args.get("fgi_idx", 0))

    def cmd_try_precursors(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """沙盒试 LLM 自提前体。args: precursors (list), reaction_type"""
        precursors = args.get("precursors", [])
        if not precursors:
            return {"error": "precursors list required"}
        return self.session.sandbox_try(
            precursors=precursors,
            reaction_type=args.get("reaction_type", ""),
        )

    def cmd_sandbox_list(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """查看沙盒历史。"""
        attempts = self.session._sandbox_attempts
        selected = self.session._sandbox_selected
        return {
            "n_attempts": len(attempts),
            "selected": selected,
            "attempts": [
                {
                    "idx": i,
                    "source": a.get("source", "?"),
                    "success": a.get("success", False),
                    "precursors": a.get("precursors", []),
                    "reaction_type": a.get("reaction_type", ""),
                    "forward_pass": (a.get("forward_validation") or {}).get("pass"),
                }
                for i, a in enumerate(attempts)
            ],
        }

    def cmd_sandbox_clear(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """清空沙盒。"""
        self.session.sandbox_clear()
        return {"ok": True, "hint": "沙盒已清空"}

    def cmd_select(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """选中沙盒方案。args: idx"""
        idx = args.get("idx", 0)
        return self.session.sandbox_select(idx)

    # ── 决策提交 ──

    def cmd_commit(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """提交决策写入树。
        args: idx (attempt), reasoning, confidence, rejected (list)
        """
        return self.session.commit(
            attempt_idx=args.get("idx"),
            reasoning=args.get("reasoning", ""),
            confidence=args.get("confidence", "medium"),
            rejected=args.get("rejected"),
        )

    def cmd_accept(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """标记当前分子为 terminal。args: reason"""
        self.session.accept_terminal(reason=args.get("reason", ""))
        return {"ok": True, "action": "accepted_terminal"}

    def cmd_skip(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """跳过当前分子。args: reason"""
        self.session.skip(reason=args.get("reason", ""))
        return {"ok": True, "action": "skipped"}

    # ── 查看 ──

    def cmd_tree(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """打印合成树。"""
        tree_text = self.session.orch.tree.print_tree()
        terminals = self.session.orch.tree.get_terminal_smiles()
        pending = self.session.orch.tree.get_pending_molecules()
        return {
            "tree": tree_text,
            "terminal_count": len(terminals),
            "pending_count": len(pending),
            "terminals": terminals,
        }

    def cmd_status(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """查看编排状态。"""
        return self.session.orch.get_status()

    def cmd_finalize(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """完成编排。args: summary"""
        return self.session.finalize(summary=args.get("summary", ""))

    def cmd_report(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """生成正向合成报告。"""
        tree = self.session.orch.tree
        report_text = generate_forward_report(tree)
        terminals = get_terminal_list(tree)
        return {
            "report": report_text,
            "starting_materials": terminals,
        }

    def cmd_export(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """导出结果到 output/ 目录。args: name, output_dir"""
        return export_results(
            self.session,
            output_dir=args.get("output_dir"),
            name=args.get("name"),
        )

    def cmd_smart_cap(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """智能断键推理。args: bond_idx 或 bond (atom pair) + smiles (可选)"""
        from Rachel.chem_tools.smart_cap import suggest_capping

        # 方式 1: 直接指定 smiles + bond
        smiles = args.get("smiles")
        bond = args.get("bond")
        if smiles and bond:
            return suggest_capping(smiles, tuple(bond), max_proposals=args.get("max", 5))

        # 方式 2: 从当前 context 中按 bond_idx 查询
        bond_idx = args.get("bond_idx", 0)
        orch = self.session.orch
        ctx = orch._current_context
        if ctx is None or ctx.decision_context is None:
            return {"error": "no active context, use smiles+bond or call next first"}

        bonds = ctx.decision_context.get("disconnectable_bonds", [])
        if bond_idx < 0 or bond_idx >= len(bonds):
            return {"error": "bond_idx %d out of range (0-%d)" % (bond_idx, len(bonds) - 1)}

        atoms = tuple(bonds[bond_idx]["atoms"])
        return suggest_capping(ctx.smiles, atoms, max_proposals=args.get("max", 5))

    def cmd_custom_cap(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """LLM 自定义 capping。args: cap_i, cap_j + (bond_idx 或 smiles+bond), reaction_type"""
        from Rachel.chem_tools.smart_cap import custom_cap

        cap_i = args.get("cap_i", "")
        cap_j = args.get("cap_j", "")
        if not cap_i or not cap_j:
            return {"error": "cap_i and cap_j required"}

        # 方式 1: 直接指定 smiles + bond
        smiles = args.get("smiles")
        bond = args.get("bond")
        if smiles and bond:
            return custom_cap(smiles, tuple(bond), cap_i, cap_j,
                              reaction_type=args.get("reaction_type", "llm_custom"))

        # 方式 2: 从当前 context 按 bond_idx
        bond_idx = args.get("bond_idx", 0)
        orch = self.session.orch
        ctx = orch._current_context
        if ctx is None or ctx.decision_context is None:
            return {"error": "no active context, use smiles+bond or call next first"}

        bonds = ctx.decision_context.get("disconnectable_bonds", [])
        if bond_idx < 0 or bond_idx >= len(bonds):
            return {"error": "bond_idx %d out of range (0-%d)" % (bond_idx, len(bonds) - 1)}

        atoms = tuple(bonds[bond_idx]["atoms"])
        return custom_cap(ctx.smiles, atoms, cap_i, cap_j,
                          reaction_type=args.get("reaction_type", "llm_custom"))


# ── 固定路径 JSON 通信 ──
# .rachel/ 目录下:
#   cmd.json     — LLM 写入命令 {"command": "next", "args": {...}}
#   result.json  — 系统写入结果
#   session.json — 会话状态（持久化）
#
# 两种 CLI 模式:
#   模式 1 (直接): python retro_cmd.py <command> [session] [args.json]
#   模式 2 (JSON): python retro_cmd.py --run [rachel_dir]
#     读 .rachel/cmd.json → 执行 → 写 .rachel/result.json → 清空 cmd.json

RACHEL_DIR = Path(__file__).resolve().parent.parent / ".rachel"
DEFAULT_SESSION = str(RACHEL_DIR / "session.json")
CMD_FILE = RACHEL_DIR / "cmd.json"
RESULT_FILE = RACHEL_DIR / "result.json"


def _ensure_rachel_dir():
    """确保 .rachel/ 目录存在。"""
    RACHEL_DIR.mkdir(parents=True, exist_ok=True)


def _run_from_cmd_json(rachel_dir: Optional[Path] = None):
    """模式 2: 读 cmd.json → 执行 → 写 result.json。"""
    rdir = rachel_dir or RACHEL_DIR
    cmd_file = rdir / "cmd.json"
    result_file = rdir / "result.json"
    session_file = str(rdir / "session.json")

    if not cmd_file.exists():
        result = {"error": f"cmd.json not found at {cmd_file}"}
        result_file.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
        return

    with open(cmd_file, "r", encoding="utf-8") as f:
        cmd_data = json.load(f)

    command = cmd_data.get("command", "")
    args = cmd_data.get("args", {})

    if not command:
        result = {"error": "cmd.json missing 'command' field"}
        result_file.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
        return

    cmd = RetroCmd(session_file)
    result = cmd.execute(command, args)

    result_file.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")

    # 清空 cmd.json（标记已消费）
    cmd_file.write_text("{}", encoding="utf-8")


def _cli_main():
    import sys
    args_list = sys.argv[1:]

    # 静默 RDKit 警告
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        pass

    _ensure_rachel_dir()

    if not args_list:
        print(json.dumps({
            "error": "usage: retro_cmd <command> [session_file] [args.json]",
            "alt_usage": "retro_cmd --run [rachel_dir]  (读 cmd.json 模式)",
            "commands": [
                "init", "next", "context", "explore", "explore_fgi",
                "try_bond", "try_fgi", "try_precursors",
                "sandbox_list", "sandbox_clear", "select",
                "commit", "accept", "skip",
                "tree", "status", "finalize", "report", "export",
                "smart_cap",
                "custom_cap",
            ],
            "rachel_dir": str(RACHEL_DIR),
            "session_file": DEFAULT_SESSION,
        }, indent=2, ensure_ascii=False))
        return

    # 模式 2: --run
    if args_list[0] == "--run":
        rdir = Path(args_list[1]) if len(args_list) > 1 else None
        _run_from_cmd_json(rdir)
        return

    # 模式 1: 直接命令
    command = args_list[0]
    session_file = args_list[1] if len(args_list) > 1 else DEFAULT_SESSION
    args_input = args_list[2] if len(args_list) > 2 else None

    # 解析 args
    cmd_args: Dict[str, Any] = {}
    if args_input:
        args_path = Path(args_input)
        if args_path.exists() and args_path.suffix == ".json":
            with open(args_path, "r", encoding="utf-8") as f:
                cmd_args = json.load(f)
        else:
            try:
                cmd_args = json.loads(args_input)
            except json.JSONDecodeError:
                print(json.dumps({"error": f"invalid args: {args_input}"}))
                return

    cmd = RetroCmd(session_file)
    result = cmd.execute(command, cmd_args)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    _cli_main()
