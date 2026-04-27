# 1word-2^61-1

固定 modulus `B = 2^61 - 1` (Mersenne 素数) に対し、`a % B` を `N` 個の `a` について計算するベンチ問題。`B` はコンパイル時定数。

各 `a_i` は base.cpp 内で `(sa, ma, mb)` から決定的に生成する:

```
sa <- sa * ma + mb  (u64 wrap)
a  <- sa & ((1 << 63) - 1)
acc ^= div.mod(a)
```

`B` がコンパイル時に既知なので、コンパイラ自身が naive `a % B` を multiply-by-magic に置き換える。Mersenne 素数特有の shift-add トリック (`(a & B) + (a >> 61)`) と、汎用 Barrett / Anton 等を比較する。

`a < 2^63` を保証しているため `(a & B) + (a >> 61)` は高々 1 回の補正で `[0, B)` に収まる。

## 入力

1 行に 4 整数を与える。

`N sa ma mb`

- `0 <= N`
- `0 <= sa, ma, mb < 2^64`

## 出力

`acc` (各 `a_i` に対する `B` での剰余の累積 XOR) を出力する。

## ケース方針

- 小さい手書きケース
- 様々な N でのランダムケース

## 再生成

以下を問題ディレクトリで実行する。

```bash
./gen/generate.sh
```
