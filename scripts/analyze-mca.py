#!/usr/bin/env python3
# /// script
# requires-python = ">=3.12"
# ///
"""llvm-mca による hot function の静的 ILP 解析。

Usage:
  python3 analyze-mca.py <binary> <callgrind_out> [--mcpu MCPU]

入力:
  - binary:        対象 ELF/Mach-O
  - callgrind_out: callgrind --dump-instr=yes 出力 (run-callgrind.sh と同設定)

処理:
  1. callgrind の per-function Ir 集計から **最も Ir の多い関数**を選ぶ
     (libc/標準ライブラリ等のシンボルは除外)
  2. objdump --disassemble=<sym> で当該関数のアセンブリを取得
  3. アドレス/オペコードバイトを除去して llvm-mca に流す
  4. JSON 出力 (--json) から IPC / Block RThroughput / Bottlenecks を抽出

出力 (stdout に 1 行 JSON):
  {
    "mca_function": "yosupo_283661::isprime_impl::isprime",
    "mca_ipc": 3.42,
    "mca_block_rthroughput": 12.5,
    "mca_uops_per_cycle": 3.50,
    "mca_bottleneck": "Resource pressure: Port 1, Port 5",
    "mca_mcpu": "skylake-avx512"
  }

llvm-mca が PATH にない / 解析失敗の場合は静的に判明する fields だけ出す。
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# 「ユーザコード」とみなす shared object 名のヒント。fn= の前に ob= で出る。
# 全部の libc 系を除外したい。判定は単に「ob= が binary 自身でなければ skip」。
_LIB_OB_RE = re.compile(r"(libc|libstdc\+\+|ld-linux|libpthread|libm|libgcc|libdyld|libsystem)")


def _find_llvm_mca() -> str | None:
    """llvm-mca バイナリのパスを探す。

    Ubuntu の `apt install llvm` は `/usr/bin/llvm-mca-NN` (versioned) しか
    入らないことがあるので、unversioned で見つからなければ versioned も探す。
    """
    p = shutil.which("llvm-mca")
    if p:
        return p
    # Ubuntu noble は llvm-18 がデフォルト。20 → 14 まで降順に探索。
    for v in range(20, 13, -1):
        p = shutil.which(f"llvm-mca-{v}")
        if p:
            return p
    # Homebrew (macOS)
    for cand in ("/opt/homebrew/opt/llvm/bin/llvm-mca",
                 "/usr/local/opt/llvm/bin/llvm-mca"):
        if Path(cand).exists():
            return cand
    return None


def _detect_arch(binary: str) -> str:
    try:
        out = subprocess.run(["file", "-b", binary],
                             capture_output=True, text=True, check=True).stdout.lower()
    except Exception:
        return "unknown"
    if "x86-64" in out or "x86_64" in out:
        return "x86_64"
    if "aarch64" in out or "arm64" in out:
        return "aarch64"
    return "unknown"


def _default_mcpu(arch: str) -> str:
    """GHA ランナー / 開発機の代表 CPU。実 CPU と完全一致しなくても傾向は出る。"""
    if arch == "x86_64":
        return "skylake-avx512"  # GHA ubuntu-latest は Xeon Cascade Lake 系
    if arch == "aarch64":
        return "neoverse-n1"  # GHA arm runner / Graviton 系
    return "generic"


_HEADER_PREFIXES = (
    "fn=", "fl=", "ob=", "cfn=", "cfi=", "cob=", "cfl=", "calls=",
    "positions:", "events:", "version:", "creator:", "desc:",
    "summary:", "totals:", "pid:", "cmd:", "part:", "thread:", "#",
)

# callgrind の "fn=" 値: 初回は `fn=(123) actual_name` (compressed-functions=yes デフォルト)。
# 2 回目以降は `fn=(123)` だけ (= name は省略)。先頭の `(NNNN) ` を除いた本体を抽出する。
# 本体が空なら ID のみで、グローバルなテーブルから引き直す必要がある。
_FN_ID_RE = re.compile(r"^\((\d+)\)\s*(.*)$")

# 実関数として扱わない hot function 候補。これらは callgrind の attribution 集約点で、
# 我々が逆アセンブルして MCA したい対象ではない。
# 注意: `main` は除外しない。-O2 で最適化された binary では、ユーザコード関数が
# 全部 main に inline 化されるケースがあり、その場合 main を逆アセンブルするのが正解。
_FN_EXCLUDE = {
    "(below main)", "(below ???)", "???",
    "_start", "_init", "_fini", "__libc_start_main",
    "__libc_csu_init", "__libc_csu_fini", "_dl_relocate_static_pie",
}
_FN_EXCLUDE_PREFIXES = ()  # 個別名のみ除外、prefix 除外は使わない


def _resolve_fn_value(value: str, fn_table: dict[str, str]) -> str:
    """callgrind の fn= 値を実名に解決する。

    `(NNNN) name` 形式なら name を fn_table[NNNN] に登録して name を返す。
    `(NNNN)` だけなら fn_table から引いて返す。
    数字 prefix が無いなら value をそのまま返す。
    """
    m = _FN_ID_RE.match(value)
    if not m:
        return value
    fid, name = m.group(1), m.group(2).strip()
    if name:
        fn_table[fid] = name
        return name
    # ID のみ: 既出のはずだから fn_table から引く。失敗したら value そのまま。
    return fn_table.get(fid, value)


def _is_excluded_fn(name: str) -> bool:
    if name in _FN_EXCLUDE:
        return True
    # main / main'2 / main.cold / __libc_start_main 等。'(...)' 全体が除外名と一致するもの。
    base = name.split("'", 1)[0].split(".", 1)[0]
    if base in _FN_EXCLUDE_PREFIXES:
        return True
    return False


def _parse_hot_function(cg_path: str, binary_path: str) -> str | None:
    """callgrind 出力から「ユーザコード関数で Ir 最大」を返す。

    バイナリ自身の ob= 行に属する fn= ブロックのみを対象、かつ
    callgrind の attribution 集約用エントリ ((below main) 等) は除外する。
    """
    binary_basename = Path(binary_path).name
    fn_ir: dict[str, int] = defaultdict(int)
    fn_table: dict[str, str] = {}  # callgrind compressed function id → 実名
    cur_ob: str = ""
    cur_fn: str = ""
    events: list[str] = []
    positions: list[str] = ["line"]
    with open(cg_path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith("ob="):
                cur_ob = line[3:].strip()
                continue
            if line.startswith("fn="):
                cur_fn = _resolve_fn_value(line[3:].strip(), fn_table)
                continue
            # cfn=(NNNN) name 形式でも name 登録が起きるので、テーブルを更新する。
            # cfn= 自体は call relationship なので Ir 集計には使わない (fn= とは別)。
            if line.startswith("cfn="):
                _resolve_fn_value(line[4:].strip(), fn_table)
                continue
            if line.startswith("positions:"):
                positions = line[len("positions:"):].split()
                continue
            if line.startswith("events:"):
                events = line[len("events:"):].split()
                continue
            if any(line.startswith(p) for p in _HEADER_PREFIXES):
                continue
            if not cur_fn or not cur_ob:
                continue
            if binary_basename not in cur_ob:
                continue
            if _LIB_OB_RE.search(cur_ob):
                continue
            if _is_excluded_fn(cur_fn):
                continue
            try:
                ir_idx = events.index("Ir")
            except ValueError:
                continue
            parts = line.split()
            try:
                ir = int(parts[len(positions) + ir_idx])
            except (ValueError, IndexError):
                continue
            fn_ir[cur_fn] += ir
    if not fn_ir:
        return None
    return max(fn_ir.items(), key=lambda kv: kv[1])[0]


def _disassemble_function(binary: str, fn: str, arch: str) -> str | None:
    """objdump でデマングル + 該当関数のアセンブリを返す (header / address 除去済み)。"""
    # mangling されたシンボル名が必要なので、まずシンボル一覧から fn にマッチするものを探す。
    # objdump --syms で探すより、デマングル後の名前で全て眺めて候補を採る。
    try:
        nm = subprocess.run(
            ["objdump", "--syms", "-C", binary],
            capture_output=True, text=True, check=True,
        ).stdout
    except Exception as e:
        print(f"objdump --syms failed: {e}", file=sys.stderr)
        return None
    # 一致条件: シンボル名に fn が含まれる (デマングル後)
    candidates: list[str] = []
    for line in nm.splitlines():
        if fn in line:
            # 行末がシンボル名。空白区切りで最後のフィールド (デマングル名は空白を含むので)
            # ここではシンプルに「fn を含む行」を全部拾う。
            # マングルされた raw 名は --syms の前半に出るので、それを別途引く。
            pass
    # mangled name を引くために --syms (no -C)
    try:
        nm_raw = subprocess.run(
            ["objdump", "--syms", binary],
            capture_output=True, text=True, check=True,
        ).stdout
    except Exception as e:
        print(f"objdump --syms (raw) failed: {e}", file=sys.stderr)
        return None
    raw_to_demangled: dict[str, str] = {}
    for raw_line, dem_line in zip(nm_raw.splitlines(), nm.splitlines()):
        toks_raw = raw_line.split()
        toks_dem = dem_line.split(maxsplit=len(raw_line.split()) - 1)
        if not toks_raw:
            continue
        raw_sym = toks_raw[-1]
        # demangled name は同じ位置 (last field の代わりに残り全部)
        if len(toks_dem) >= len(toks_raw):
            dem_sym = toks_dem[-1]
        else:
            dem_sym = raw_sym
        raw_to_demangled[raw_sym] = dem_sym
    # fn を含む raw シンボルを採用 (マングル名 OR デマングル名どちらかで一致)
    for raw, dem in raw_to_demangled.items():
        if fn == raw or fn == dem or fn in dem:
            candidates.append(raw)
    if not candidates:
        # callgrind が報告する fn 名と objdump --syms の出力ズレに備え、緩いマッチ
        for raw, dem in raw_to_demangled.items():
            if any(part and part in dem for part in fn.split("::")):
                candidates.append(raw)
                break
    if not candidates:
        print(f"no symbol candidate for fn={fn!r}", file=sys.stderr)
        return None
    sym = candidates[0]
    # GNU binutils (Linux) と Apple LLVM (macOS) でフラグ名が違う:
    #   GNU:   --disassemble=<sym>
    #   Apple: --disassemble-symbols=<sym>
    # 順に試す。
    dump: str | None = None
    for flag in (f"--disassemble={sym}", f"--disassemble-symbols={sym}"):
        try:
            r = subprocess.run(
                ["objdump", "-d", "--no-show-raw-insn", flag, binary],
                capture_output=True, text=True,
            )
            if r.returncode == 0 and r.stdout.strip():
                dump = r.stdout
                break
        except Exception:
            continue
    if dump is None:
        print(f"objdump disassemble failed for {sym}", file=sys.stderr)
        return None
    # アドレスを除去して命令だけ取り出す
    out_lines: list[str] = []
    started = False
    for line in dump.splitlines():
        # ヘッダ・空行スキップ
        s = line.strip()
        if not s:
            continue
        if s.startswith("Disassembly") or s.endswith(":"):
            # "<sym>:" や "Disassembly of section ..."
            started = s.endswith(":")
            continue
        if not started:
            continue
        # "  401234: mnemonic operands"  →  "mnemonic operands"
        m = re.match(r"^\s*[0-9a-fA-F]+:\s+(.*)$", line)
        if m:
            instr = m.group(1).rstrip()
            if instr and not instr.startswith("("):
                out_lines.append(instr)
    if not out_lines:
        return None
    return "\n".join(out_lines) + "\n"


def _run_llvm_mca(asm: str, mcpu: str, arch: str, mca_path: str) -> dict | None:
    """llvm-mca を asm に対して実行し、JSON を返す。"""
    triple = {
        "x86_64": "x86_64-linux-gnu",
        "aarch64": "aarch64-linux-gnu",
    }.get(arch, "")
    cmd = [mca_path, f"-mcpu={mcpu}", "--json"]
    if triple:
        cmd.append(f"-mtriple={triple}")
    try:
        r = subprocess.run(cmd, input=asm, capture_output=True, text=True, timeout=60)
    except subprocess.TimeoutExpired:
        print("llvm-mca timeout", file=sys.stderr)
        return None
    if r.returncode != 0:
        print(f"llvm-mca exit={r.returncode}\nstderr:\n{r.stderr[:500]}", file=sys.stderr)
        return None
    try:
        return json.loads(r.stdout)
    except json.JSONDecodeError as e:
        print(f"llvm-mca JSON parse error: {e}", file=sys.stderr)
        return None


def _extract_metrics(mca_json: dict) -> dict:
    """llvm-mca --json から 主要メトリクスを抽出。"""
    out: dict = {}
    # llvm-mca の JSON は "CodeRegions" -> [ { "SummaryView": {...}, ... } ] の形
    regions = mca_json.get("CodeRegions") or []
    if not regions:
        return out
    summary = regions[0].get("SummaryView") or {}
    if "IPC" in summary:
        out["mca_ipc"] = float(summary["IPC"])
    if "BlockRThroughput" in summary:
        out["mca_block_rthroughput"] = float(summary["BlockRThroughput"])
    if "uOpsPerCycle" in summary:
        out["mca_uops_per_cycle"] = float(summary["uOpsPerCycle"])
    if "Iterations" in summary:
        out["mca_iterations"] = int(summary["Iterations"])
    if "Instructions" in summary:
        out["mca_instructions"] = int(summary["Instructions"])
    # bottleneck info があれば
    bn = regions[0].get("BottleneckAnalysis") or {}
    if isinstance(bn, dict):
        # llvm-mca の BottleneckAnalysis の構造はバージョン依存だが ResourcePressure 系を文字列化
        keys = [k for k, v in bn.items() if isinstance(v, (int, float)) and v > 0]
        if keys:
            out["mca_bottleneck_keys"] = keys[:5]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("callgrind_out")
    ap.add_argument("--mcpu", default=None)
    ap.add_argument("--dry-run", action="store_true",
                    help="llvm-mca を呼ばずに hot function 抽出までで止める")
    args = ap.parse_args()

    arch = _detect_arch(args.binary)
    mcpu = args.mcpu or _default_mcpu(arch)

    fn = _parse_hot_function(args.callgrind_out, args.binary)
    out: dict = {"mca_arch": arch, "mca_mcpu": mcpu}
    if fn:
        out["mca_function"] = fn
    else:
        print(json.dumps(out))
        return

    asm = _disassemble_function(args.binary, fn, arch)
    if not asm:
        print(json.dumps(out))
        return
    out["mca_function_size_lines"] = asm.count("\n")

    if args.dry_run:
        print(json.dumps(out))
        return

    mca_path = _find_llvm_mca()
    if not mca_path:
        out["mca_error"] = "llvm-mca not found in PATH"
        print(json.dumps(out))
        return

    mca_json = _run_llvm_mca(asm, mcpu, arch, mca_path)
    if mca_json is None:
        out["mca_error"] = "llvm-mca failed (see stderr)"
        print(json.dumps(out))
        return

    out.update(_extract_metrics(mca_json))
    print(json.dumps(out))


if __name__ == "__main__":
    main()
