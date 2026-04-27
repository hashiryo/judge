# judge

競プロの作業用ジャッジシステム。同じ問題に対する複数の algo 実装を共通の計測ハーネスで比較できる。

## ディレクトリ構成

```
problems/
  yosupo-unionfind/
    problem.toml       # 問題 URL や TLE 等の設定
    base.cpp           # 計測ハーネス (入出力 + 計測ループ)
    algos/
      _common.hpp      # algo 共通ヘッダ (typedef 等。_ 始まりは提出対象外)
      naive.hpp        # algo 実装。base.cpp が ALGO_HPP 経由で include する
      rank.hpp
      simd.hpp
    testcases/         # サンプルケース (fetch-samples.py で取得、コミット対象)
      sample_00.in
      sample_00.out
.results/
  history.jsonl        # 全履歴 (GHA が自動更新)
lib/
  include/
    bits/stdc++.h      # ローカルテスト用 (macOS clang++ 向け)
    simde/             # ARM 上で x86 SIMD intrinsics をエミュレート
scripts/
  fetch-samples.py     # サンプルケースのダウンロード
  local-test.sh        # ローカルでのサンプルテスト
  detect-changed.py    # GHA: 変更された問題の検出
  download-testcases.py # GHA: 全テストケースのダウンロード
  run-tests.sh         # GHA: テスト実行
  run-callgrind.sh     # GHA: callgrind による命令数計測
  collect-results.py   # GHA: 結果を history.jsonl に追記
  write-summary.py     # GHA: Step Summary 生成
```

## ハーネス形式 (新形式)

各問題は **`base.cpp` + `algos/*.hpp`** で構成する。

- **`base.cpp`** は計測フレーム。入出力・計測ループ・`ALGO_TIME_NS` の出力を担当。
  プリプロセッサで `#include ALGO_HPP` する形になっており、CI 側で
  `-DALGO_HPP="\"algos/xxx.hpp\""` を渡すことで algo を切り替える。
- **`algos/*.hpp`** が提出候補 (= 比較対象の algo 実装)。各ファイルは `struct MP` 等の
  共通インターフェースを実装する。
- **`algos/_*.hpp`** は提出対象から除外される (`_common.hpp` 等の共通ヘッダ用)。
- 1 つの問題内で同じハーネスを使うので、入出力や計測タイミングが揃い、
  algo 間の純粋な比較ができる。

例: [problems/yosupo-unionfind/base.cpp](problems/yosupo-unionfind/base.cpp) と
[problems/yosupo-unionfind/algos/](problems/yosupo-unionfind/algos/)。

### `algo_time_ns` 出力

`base.cpp` は計測対象ループの実行時間を `chrono::steady_clock` で測り、
`stderr` に `ALGO_TIME_NS=<ナノ秒>` の行を出力する。これがケース毎に集計され、
ダッシュボードでは `algo_time_max` / `algo_time_total` として比較できる。

## セットアップ

サブモジュール (Library) とその入れ子 (simde) を再帰的に取得する:

```bash
git submodule update --init --recursive
```

ローカルテストには `timeout` コマンドが必要。macOS は Homebrew の coreutils を入れる
(入っていれば `gtimeout` として検出される):

```bash
brew install coreutils
```

## 使い方

### 1. 問題ディレクトリを作る

`problems/` 以下に問題ディレクトリを作り、`problem.toml` と `base.cpp` と `algos/*.hpp` を置く。

問題ディレクトリは `problem.toml` を持つディレクトリとして判定される。
入れ子にしてもよいが、その場合も各問題ルートには `problem.toml` が必要。

```toml
# problem.toml
url = "https://judge.yosupo.jp/problem/unionfind"
tle = 5
# error = 1e-6                  # 浮動小数点許容誤差
# callgrind_case = "mid_00"     # callgrind 用代表ケース (デフォルトは中央サイズ)
```

### 2. サンプルケースを取得

```bash
python3 scripts/fetch-samples.py problems/yosupo-unionfind/
```

yosupo / AOJ / yukicoder に対応。取得したサンプルは `testcases/` に保存されるのでコミットしておく。

### 3. ローカルで確認

```bash
./scripts/local-test.sh problems/yosupo-unionfind/algos/naive.hpp  # 特定 algo
./scripts/local-test.sh problems/yosupo-unionfind/                  # 全 algo
```

### 4. push → GHA で本番テスト

push すると GHA が自動で全テストケースをダウンロードし、4 環境 (x64/arm × g++/clang++) で実行。
同じ問題の全提出を比較した結果が Step Summary とダッシュボード (`docs/index.html`) に表示される。

## problem.toml

各問題ディレクトリに `problem.toml` を置いて設定する:

```toml
url = "https://judge.yosupo.jp/problem/unionfind"  # 問題 URL
tle = 5                                              # TLE 秒数 (デフォルト: 10)
mle = 256                                            # MLE メガバイト (デフォルト: 256)
# error = 1e-6                                       # 浮動小数点許容誤差
# callgrind_case = "mid_00"                          # callgrind 用ケース指定
```

## 自作テストケース

`testcases/` に `.in` / `.out` ペアを置けば、外部ジャッジの問題と併用できる:

```
problems/my-problem/
  problem.toml
  base.cpp
  algos/
    _common.hpp
    naive.hpp
  testcases/
    00.in  00.out
    01.in  01.out
```

カテゴリ分けのために入れ子にしてもよい:

```
problems/graph/unionfind/
  problem.toml
  base.cpp
  algos/
    naive.hpp
  testcases/
    00.in  00.out
```

`gen/` 以下はテストケース生成器の置き場として、提出候補の検出からも除外される。

## 履歴比較とアーカイブ

ダッシュボードの History タブはコミット間の `algo_time_max` の推移を表示する。
比較は **同じ問題定義 (`cases_hash`)** の範囲内に閉じられる:

- `cases_hash` は `testcases/*.in` の内容 + `base.cpp` + `problem.toml` から計算
- これらが変わると "問題定義の変更" とみなされ、それ以前の履歴は archived 扱いになる
- `algos/*.hpp` (algo 実装) と `algos/_*.hpp` (algo 共通ヘッダ) は `cases_hash` に含めない
  (algo 実装の改善やリグレッション比較を妨げないため)
- マシンタイプ (`cpu_model`) ごとにグループを分けて表示するので、
  GHA ランナーの差異が比較を歪めないようになっている
