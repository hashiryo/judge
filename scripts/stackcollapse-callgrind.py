#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""Callgrind の出力を FlameGraph 用 folded stack 形式に変換する。

Usage:
  python3 stackcollapse-callgrind.py callgrind.out > folded.txt
  perl flamegraph.pl folded.txt > flame.svg

Callgrind format (callgrind -> self/inclusive Ir):
  fn=(id) [name]        関数ブロック開始
  <line> <ir>           現在関数の self Ir
  cfn=(id) [name]       被呼関数
  calls=<n> <line>      calls 行 (直後の cost 行が inclusive Ir)
  <line> <ir>           直前の calls に対する inclusive Ir

各関数の self Ir を、incoming call の inclusive Ir 比で parent へ配分して折り畳む。
再帰はノード seen セットで打ち切り。
"""
from __future__ import annotations
import re
import sys
from collections import defaultdict


_FN_RE = re.compile(r"^(fn|cfn)=\((\d+)\)(?:\s+(.+))?$")
_OB_RE = re.compile(r"^(ob|cob)=\((\d+)\)(?:\s+(.+))?$")
_FL_RE = re.compile(r"^(fl|cfi|cfl)=\((\d+)\)(?:\s+(.+))?$")


def parse_callgrind(path: str):
    fn_names: dict[int, str] = {}
    fn_self: dict[int, int] = defaultdict(int)
    # fn_calls[caller][callee] = inclusive Ir (集計済み)
    fn_calls: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))

    current_fn: int | None = None
    pending_cfn: int | None = None
    next_cost_is_inclusive = False

    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("desc:") or line.startswith("summary:") or line.startswith("totals:") or line.startswith("events:") or line.startswith("positions:") or line.startswith("version:") or line.startswith("creator:") or line.startswith("pid:") or line.startswith("cmd:") or line.startswith("part:"):
                continue

            # fn=(id) [name] / cfn=
            m = _FN_RE.match(line)
            if m:
                kind, fid_s, name = m.group(1), m.group(2), m.group(3)
                fid = int(fid_s)
                if name:
                    fn_names[fid] = name
                if kind == "fn":
                    current_fn = fid
                    pending_cfn = None
                    next_cost_is_inclusive = False
                else:  # cfn
                    pending_cfn = fid
                continue

            # ob / cob / fl / cfi / cfl are context; ignore names but skip parse errors
            if _OB_RE.match(line) or _FL_RE.match(line):
                continue

            if line.startswith("calls="):
                next_cost_is_inclusive = True
                continue

            # cost line: "<line> <ir>" (またはそれ以上の列、末尾が cost)
            parts = line.split()
            if not parts:
                continue
            # 末尾を Ir とみなす
            last = parts[-1]
            if not (last.lstrip("+-").isdigit()):
                # フォーマット外の行は無視
                continue
            ir = int(last)

            if next_cost_is_inclusive:
                if current_fn is not None and pending_cfn is not None:
                    fn_calls[current_fn][pending_cfn] += ir
                next_cost_is_inclusive = False
                pending_cfn = None
            else:
                if current_fn is not None:
                    fn_self[current_fn] += ir

    return fn_names, dict(fn_self), {k: dict(v) for k, v in fn_calls.items()}


def fold(fn_names, fn_self, fn_calls):
    # callers[callee][caller] = inclusive Ir
    callers: dict[int, dict[int, int]] = defaultdict(dict)
    for caller, callees in fn_calls.items():
        for callee, ir in callees.items():
            callers[callee][caller] = ir

    def label(fn: int) -> str:
        name = fn_names.get(fn, f"fn_{fn}")
        # FlameGraph folded format separator を避ける
        return name.replace(";", ":").replace(" ", "_")

    # root 判定: 呼び出し元が無い関数 (self-recursive も root にする)
    all_fns = set(fn_self.keys()) | set(fn_calls.keys()) | set(callers.keys())

    # path の weight を流す DFS。seen で cycle を打ち切る。
    # paths_for[fn] = [(path_of_ids, weight), ...]
    memo: dict[int, list[tuple[tuple[int, ...], float]]] = {}

    def paths_for(fn: int, seen: frozenset[int]) -> list[tuple[tuple[int, ...], float]]:
        if fn in seen:
            return []
        # memo は seen 非依存版のみ (頻出パターンへの最適化)
        if not seen and fn in memo:
            return memo[fn]
        parents = callers.get(fn, {})
        if not parents:
            result = [((fn,), 1.0)]
        else:
            total = sum(parents.values())
            result = []
            new_seen = seen | {fn}
            for p, pir in parents.items():
                w = (pir / total) if total > 0 else (1.0 / len(parents))
                for parent_path, parent_w in paths_for(p, new_seen):
                    result.append((parent_path + (fn,), parent_w * w))
            if not result:
                # すべて cycle だった場合は self を root 扱い
                result = [((fn,), 1.0)]
        if not seen:
            memo[fn] = result
        return result

    out_lines: list[str] = []
    for fn in all_fns:
        s = fn_self.get(fn, 0)
        if s <= 0:
            continue
        paths = paths_for(fn, frozenset())
        for path, w in paths:
            contribution = int(s * w)
            if contribution <= 0:
                continue
            out_lines.append(";".join(label(x) for x in path) + " " + str(contribution))
    return out_lines


def main():
    if len(sys.argv) < 2:
        print("Usage: stackcollapse-callgrind.py <callgrind.out>", file=sys.stderr)
        sys.exit(1)
    fn_names, fn_self, fn_calls = parse_callgrind(sys.argv[1])
    for line in fold(fn_names, fn_self, fn_calls):
        print(line)


if __name__ == "__main__":
    main()
