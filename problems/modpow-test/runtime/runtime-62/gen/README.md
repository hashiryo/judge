# runtime-62 (modpow)

入力で法 `mod` を実行時に与える、64-bit 動的 mod での冪乗ベンチ問題。

各クエリ `(a_i, b_i)` に対して `pow(a_i, b_i) mod mod` を計算し、その結果の XOR を返す。

## 入力

1 行に 5 整数を与える。

`N mod bbits sa sb`

- `0 <= N`
- `2^40 < mod < 2^62`
- `mod` は奇数
- `1 <= bbits <= 64`（指数 `b_i` のビット幅）
- `0 <= sa, sb < 2^64`（LCG の初期 seed）

クエリは Knuth LCG（`s <- s * 6364136223846793005 + 1442695040888963407`）を 2 系列回し、
`a_i = sa_i % mod`, `b_i = sb_i & ((1 << bbits) - 1)` として `N` 個生成する。

## 出力

`acc = XOR_{i=0..N-1} (a_i^{b_i} mod mod)` を 1 行で出力する。

## 狙い

- 内部 mul は約 `bbits * N` 回。`bbits = 60`, `N = 5 * 10^6` で `3 * 10^8` mul と十分な信号量。
- 通常の Montgomery と Plantard でどちらが速いかを比較する。
- 将来的に「base 側を preconditioned に持つ Plantard 変種で 1 mul 削れる」かを検証する。
