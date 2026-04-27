# 1word-998244353

固定 modulus `B = 998244353` (NTT 用 30 bit 素数) に対し、`a % B` を `N` 個の `a` について計算するベンチ問題。`B` はコンパイル時定数。

各 `a_i` は base.cpp 内で `(sa, ma, mb)` から決定的に生成する:

```
sa <- sa * ma + mb  (u64 wrap)
a  <- u32(sa) & ((1 << 31) - 1)
acc ^= div.mod(a)
```

`B` がコンパイル時に既知なので、コンパイラ自身が naive `a % B` を multiply-by-magic に置き換える。本問題ではそれと、明示的に書いた Barrett / Lemire / Anton 等を比較する。

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
