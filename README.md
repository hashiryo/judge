# judge

競プロの作業用ジャッジシステム。問題単位で複数の提出を比較できる。

## ディレクトリ構成

```
problems/
  yosupo-unionfind/
    problem.toml       # 問題 URL や TLE 等の設定
    sol_a.cpp
    sol_b.cpp
    testcases/         # サンプルケース (fetch-samples.py で取得、コミット対象)
      sample_00.in
      sample_00.out
.results/
  history.jsonl        # 全履歴 (GHA が自動更新)
include/
  bits/stdc++.h        # ローカルテスト用 (macOS clang++ 向け)
scripts/
  fetch-samples.py     # サンプルケースのダウンロード
  local-test.sh        # ローカルでのサンプルテスト
  detect-changed.py    # GHA: 変更された問題の検出
  download-testcases.py # GHA: 全テストケースのダウンロード
  run-tests.sh         # GHA: テスト実行
  collect-results.py   # GHA: 結果を history.jsonl に追記
  write-summary.py     # GHA: Step Summary 生成
```

## 使い方

### 1. 問題ディレクトリを作る

`problems/` 以下にディレクトリを作り、`problem.toml` と `.cpp` ファイルを置く。

```toml
# problem.toml
url = "https://judge.yosupo.jp/problem/unionfind"
tle = 5
# error = 1e-6
```

### 2. サンプルケースを取得

```bash
python3 scripts/fetch-samples.py problems/yosupo-unionfind/
```

yosupo / AOJ / yukicoder に対応。取得したサンプルは `testcases/` に保存されるのでコミットしておく。

### 3. ローカルで確認

```bash
./scripts/local-test.sh problems/yosupo-unionfind/sol_a.cpp  # 特定ファイル
./scripts/local-test.sh problems/yosupo-unionfind/            # 全ファイル
```

### 4. push → GHA で本番テスト

push すると GHA が自動で全テストケースをダウンロードし、4環境 (x64/arm × g++/clang++) で実行。
同じ問題の全提出を比較した結果が Step Summary に表示される。

## problem.toml

各問題ディレクトリに `problem.toml` を置いて設定する:

```toml
url = "https://judge.yosupo.jp/problem/unionfind"  # 問題 URL
tle = 5                                              # TLE 秒数 (デフォルト: 10)
# error = 1e-6                                       # 浮動小数点許容誤差
```

## 自作テストケース

`testcases/` に `.in` / `.out` ペアを置けば、外部ジャッジの問題と併用できる:

```
problems/my-problem/
  sol.cpp
  testcases/
    00.in  00.out
    01.in  01.out
```
