#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""Callgrind の JSONL 結果を既存 result-*.json にマージする。

入力:
  .cache/results/result-${ENV}.json    既存のテスト結果 (entry list)
  .cache/callgrind/callgrind-${ENV}.jsonl  Callgrind 結果 (JSONL, {file,environment,callgrind_*})

処理:
  (file, environment) でキー結合し、既存 entry の top-level に
  callgrind_instructions, callgrind_case を付与する。
  両側にあるファイルだけ処理される (callgrind 側に無ければスキップ)。

出力:
  元の result-${ENV}.json を上書き。
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULT_DIR = ROOT / ".cache" / "results"
CG_DIR = ROOT / ".cache" / "callgrind"


def load_callgrind_map() -> dict[tuple[str, str], dict]:
    """callgrind-*.jsonl を読み、(file, env) → {callgrind_*} のマップを返す。"""
    mapping: dict[tuple[str, str], dict] = {}
    if not CG_DIR.exists():
        return mapping
    for jsonl in sorted(CG_DIR.glob("callgrind-*.jsonl")):
        with open(jsonl) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                data = json.loads(line)
                key = (data["file"], data["environment"])
                mapping[key] = {
                    "callgrind_instructions": data.get("callgrind_instructions"),
                    "callgrind_case": data.get("callgrind_case"),
                    "callgrind_d1_misses": data.get("callgrind_d1_misses"),
                    "callgrind_ll_misses": data.get("callgrind_ll_misses"),
                    "callgrind_branch_misses": data.get("callgrind_branch_misses"),
                }
    return mapping


def merge_into_result_files(cg_map: dict[tuple[str, str], dict]) -> int:
    """result-*.json の各 entry に Callgrind 情報を付与し、上書き保存する。"""
    merged = 0
    for rf in sorted(RESULT_DIR.glob("result-*.json")):
        with open(rf) as f:
            entries = json.load(f)
        if not isinstance(entries, list):
            continue
        changed = False
        for entry in entries:
            key = (entry.get("file", ""), entry.get("environment", ""))
            cg = cg_map.get(key)
            if cg:
                for k, v in cg.items():
                    if v is not None:
                        entry[k] = v
                changed = True
                merged += 1
        if changed:
            with open(rf, "w") as f:
                json.dump(entries, f, indent=2, ensure_ascii=False)
                f.write("\n")
    return merged


def main() -> None:
    cg_map = load_callgrind_map()
    if not cg_map:
        print("No callgrind results to merge")
        return
    count = merge_into_result_files(cg_map)
    print(f"Merged {count} callgrind entries into {len(list(RESULT_DIR.glob('result-*.json')))} result files")


if __name__ == "__main__":
    main()
