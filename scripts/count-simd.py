#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""バイナリと Callgrind dump-instr 出力から SIMD 命令の Ir を集計する。

Usage:
  python3 count-simd.py <binary> <callgrind_out>

前提:
  callgrind は --dump-instr=yes --compress-pos=no で実行されていること。

SIMD 判定:
  - x86_64: ニーモニックが 'v' で始まる (AVX/AVX2/AVX-512 は v prefix)
  - aarch64: オペランドにベクトルレジスタ (v<N>.XXX または z<N>.XXX) を含む

出力:
  "<simd_ir> <total_ir>" を stdout に 1 行で出力。
  total_ir はバイナリ自身のコードに限定 (addr が disasm に含まれるもの)。
  ライブラリコード (libc 等) は total に含まない。
"""
from __future__ import annotations
import re
import subprocess
import sys


# オペランド中のベクトルレジスタ参照 (arm64 NEON/SVE 用)
_ARM_VEC_RE = re.compile(r"\b[vz]\d+\.[0-9a-z]+")


def disassemble(binary_path: str) -> dict[int, tuple[str, str]]:
    """objdump -d でバイナリを逆アセンブルし、{addr: (mnemonic, operands)} を返す。"""
    result = subprocess.run(
        ["objdump", "-d", "--no-show-raw-insn", binary_path],
        capture_output=True, text=True, check=True,
    )
    # "  400620:	add    %rax,%rbx"
    line_re = re.compile(r"^\s+([0-9a-fA-F]+):\s+(\S+)\s*(.*)$")
    out: dict[int, tuple[str, str]] = {}
    for line in result.stdout.splitlines():
        m = line_re.match(line)
        if m:
            addr = int(m.group(1), 16)
            mnemonic = m.group(2)
            operands = m.group(3).strip()
            # "<data>" 形式の疑似命令を除外 (e.g. `(bad)`)
            if mnemonic.startswith("(") or mnemonic == ".byte":
                continue
            out[addr] = (mnemonic, operands)
    return out


def is_simd(mnemonic: str, operands: str, arch: str) -> bool:
    if arch.startswith("x86") or arch == "amd64":
        # AVX 系は v 始まり。SSE も -march=x86-64-v3 では v 化される
        return mnemonic.startswith("v")
    # arm64: ベクトルレジスタを含むか
    return bool(_ARM_VEC_RE.search(operands))


def detect_arch(binary_path: str) -> str:
    try:
        out = subprocess.run(
            ["file", "-b", binary_path],
            capture_output=True, text=True, check=True,
        ).stdout
    except Exception:
        return "unknown"
    low = out.lower()
    if "x86-64" in low or "x86_64" in low:
        return "x86_64"
    if "aarch64" in low or "arm64" in low:
        return "aarch64"
    return "unknown"


def parse_callgrind_instr(path: str):
    """Callgrind --dump-instr=yes --compress-pos=no 出力から (addr, Ir) を yield する。"""
    positions: list[str] = ["line"]
    events: list[str] = []
    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("positions:"):
                positions = line[len("positions:"):].strip().split()
                continue
            if line.startswith("events:"):
                events = line[len("events:"):].strip().split()
                continue
            # メタ行はスキップ
            if (
                line.startswith("fn=") or line.startswith("fl=") or line.startswith("ob=")
                or line.startswith("cfn=") or line.startswith("cfi=") or line.startswith("cob=")
                or line.startswith("cfl=") or line.startswith("calls=")
                or line.startswith("#") or line.startswith("desc:")
                or line.startswith("summary:") or line.startswith("totals:")
                or line.startswith("version:") or line.startswith("creator:")
                or line.startswith("pid:") or line.startswith("cmd:")
                or line.startswith("part:") or line.startswith("thread:")
            ):
                continue
            if "instr" not in positions or "Ir" not in events:
                continue
            parts = line.split()
            if not parts:
                continue
            instr_idx = positions.index("instr")
            ir_idx = events.index("Ir")
            ev_start = len(positions)
            try:
                instr_field = parts[instr_idx]
                # --compress-pos=no 前提なので絶対アドレスのみ
                if instr_field.startswith("0x") or all(c in "0123456789abcdefABCDEF" for c in instr_field):
                    addr = int(instr_field, 16)
                else:
                    continue
                ir_val = int(parts[ev_start + ir_idx])
                yield (addr, ir_val)
            except (ValueError, IndexError):
                continue


def main() -> None:
    if len(sys.argv) < 3:
        print("Usage: count-simd.py <binary> <callgrind.out>", file=sys.stderr)
        sys.exit(1)
    binary, cg_out = sys.argv[1], sys.argv[2]
    arch = detect_arch(binary)
    disasm = disassemble(binary)

    simd_ir = 0
    binary_ir = 0  # バイナリに含まれる命令分のみ (ライブラリは除外)
    for addr, ir in parse_callgrind_instr(cg_out):
        insn = disasm.get(addr)
        if not insn:
            continue
        binary_ir += ir
        if is_simd(insn[0], insn[1], arch):
            simd_ir += ir

    print(f"{simd_ir} {binary_ir}")


if __name__ == "__main__":
    main()
