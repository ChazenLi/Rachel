"""Host-side Rachel command runner with append-only JSON provenance logs."""

from __future__ import annotations

import argparse
import json
import re
import sys
import traceback
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from Rachel.main.retro_cmd import RetroCmd


_INPUT_FILE_PATTERN = re.compile(
    r"^(?P<number>\d{3,})_[a-z0-9_]+\.input\.json$"
)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def _next_sequence(log_dir: Path) -> int:
    numbers = []
    if log_dir.is_dir():
        for path in log_dir.glob("*.input.json"):
            match = _INPUT_FILE_PATTERN.match(path.name)
            if match:
                numbers.append(int(match.group("number")))
    return max(numbers, default=-1) + 1


def _command_slug(command: str) -> str:
    return re.sub(r"[^a-z0-9_]+", "_", command.lower()).strip("_") or "command"


def _log_stem(log_dir: Path, command: str) -> str:
    return f"{_next_sequence(log_dir):03d}_{_command_slug(command)}"


class LoggedRetroCmd:
    """Execute Rachel commands while recording host-side input/output JSON.

    Rachel's session remains the authoritative state. These files only preserve
    the command boundary used to produce or inspect that state.
    """

    def __init__(self, session_file: str, *, log_dir: Optional[str] = None):
        self.session_file = Path(session_file)
        self.log_dir = (
            Path(log_dir) if log_dir is not None else self.session_file.parent
        )
        self._cmd = RetroCmd(str(self.session_file))

    @property
    def session(self):
        return self._cmd.session

    def execute(
        self,
        command: str,
        args: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        command_args = args or {}
        request = {"command": command, "args": command_args}
        source_dir = self.log_dir
        source_stem = _log_stem(source_dir, command)
        input_path = source_dir / f"{source_stem}.input.json"
        _write_json(input_path, request)

        try:
            result = self._cmd.execute(command, command_args)
        except Exception as exc:
            result = {
                "runner_exception": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            }
            _write_json(source_dir / f"{source_stem}.output.json", result)
            return result

        destination_dir = source_dir
        destination_stem = source_stem
        if command == "review_node" and result.get("ok"):
            active_session = result.get("active_session_file")
            if active_session:
                active_path = Path(str(active_session))
                destination_dir = active_path.parent
                if destination_dir.resolve() != source_dir.resolve():
                    destination_stem = _log_stem(destination_dir, command)
                    destination_input = (
                        destination_dir / f"{destination_stem}.input.json"
                    )
                    destination_dir.mkdir(parents=True, exist_ok=True)
                    input_path.replace(destination_input)
                self.session_file = active_path
                self.log_dir = destination_dir

        _write_json(
            destination_dir / f"{destination_stem}.output.json",
            result,
        )
        return result


def _request_from_stdin() -> Dict[str, Any]:
    request = json.load(sys.stdin)
    if not isinstance(request, dict):
        raise TypeError("stdin must contain one JSON object")
    command = request.get("command")
    args = request.get("args", {})
    if not isinstance(command, str) or not command.strip():
        raise TypeError("command must be a non-empty string")
    if not isinstance(args, dict):
        raise TypeError("args must be an object")
    return {"command": command, "args": args}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Execute one Rachel JSON command and record input/output logs."
    )
    parser.add_argument("session_file", help="Rachel session.json path")
    parser.add_argument(
        "--log-dir",
        help="Optional process-log directory; defaults to the session directory",
    )
    parsed = parser.parse_args(argv)

    try:
        request = _request_from_stdin()
    except Exception as exc:
        result = {
            "runner_exception": type(exc).__name__,
            "message": str(exc),
        }
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 1

    runner = LoggedRetroCmd(parsed.session_file, log_dir=parsed.log_dir)
    result = runner.execute(request["command"], request["args"])
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if "runner_exception" in result:
        return 1
    if "error" in result:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
