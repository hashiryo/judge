"""problem.toml を読み取る共通モジュール (簡易パーサー、stdlib のみ)。

key = "string" / key = number 形式に対応。コメント (#) と空行は無視。
"""
from __future__ import annotations
import re
from pathlib import Path


_QUOTED_RE = re.compile(r'(\w+)\s*=\s*"([^"]*)"')
_BARE_RE = re.compile(r"(\w+)\s*=\s*(\S+)")


def parse_problem_toml(problem_dir: Path) -> dict[str, str]:
    """problem.toml を読み取り、{key: str_value} の dict を返す。

    数値も str で返す点に注意 (呼び出し側で必要に応じて int/float 変換)。
    ファイルが無ければ空 dict。
    """
    toml_path = problem_dir / "problem.toml"
    config: dict[str, str] = {}
    if not toml_path.exists():
        return config
    for line in toml_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = _QUOTED_RE.match(line)
        if m:
            config[m.group(1)] = m.group(2)
            continue
        m = _BARE_RE.match(line)
        if m:
            config[m.group(1)] = m.group(2)
    return config
