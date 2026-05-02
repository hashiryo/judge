#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""Callgrind の events/totals と SIMD 集計値から結果 entry JSON を組み立てる。

入力は環境変数:
  EVENTS       "Ir Dr Dw ..." (callgrind の events ヘッダ行)
  TOTALS       "123 456 789 ..." (callgrind の totals 行)
  SIMD_IR      SIMD 命令の Ir (count-simd.py の出力)
  BINARY_IR    バイナリ本体命令の Ir (count-simd.py の出力)
  REL_PATH     sol ファイルの相対パス (ROOT 基準)
  ENV_NAME     環境名 (x64-g++ 等)
  CASE_NAME    代表ケース名
  MCA_JSON     analyze-mca.py の出力ファイルパス (任意)

出力: 1 行 JSON を stdout。
"""
import json
import os


def main() -> None:
    events = os.environ["EVENTS"].split()
    totals = list(map(int, os.environ["TOTALS"].split()))
    m = dict(zip(events, totals))
    simd_ir = int(os.environ.get("SIMD_IR", "0") or "0")
    binary_ir = int(os.environ.get("BINARY_IR", "0") or "0")
    simd_ratio = (simd_ir / binary_ir) if binary_ir > 0 else 0.0

    entry = {
        "file": os.environ["REL_PATH"],
        "environment": os.environ["ENV_NAME"],
        "callgrind_case": os.environ["CASE_NAME"],
        "callgrind_instructions": m.get("Ir", 0),
        "callgrind_d1_misses": m.get("D1mr", 0) + m.get("D1mw", 0),
        "callgrind_ll_misses": m.get("DLmr", 0) + m.get("DLmw", 0),
        "callgrind_branch_misses": m.get("Bcm", 0),
        "callgrind_simd_instructions": simd_ir,
        "callgrind_binary_instructions": binary_ir,
        "callgrind_simd_ratio": simd_ratio,
    }
    # MCA: hot function の静的 IPC 等。失敗時 (llvm-mca 不在等) は欠損で OK。
    mca_path = os.environ.get("MCA_JSON", "")
    if mca_path and os.path.exists(mca_path):
        try:
            with open(mca_path) as f:
                mca = json.load(f)
            if isinstance(mca, dict):
                # mca_* prefix の field をそのまま entry にマージ
                for k, v in mca.items():
                    if k.startswith("mca_"):
                        entry[k] = v
        except (OSError, json.JSONDecodeError):
            pass
    print(json.dumps(entry))


if __name__ == "__main__":
    main()
