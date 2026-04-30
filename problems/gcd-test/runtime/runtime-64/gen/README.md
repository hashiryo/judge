# runtime-64 (gcd)

64-bit 整数の gcd クエリベンチ問題。

各クエリ `(a_i, b_i)` に対し `gcd(a_i, b_i)` を計算し、その XOR を返す。
クエリは独立 (前のクエリの出力を入力に混ぜない) なので、OoO 実行で
algo 内のスループットが素直に出る形。

## 入力

`N abits sa sb`

- `0 <= N`
- `1 <= abits <= 64` (各 `a_i, b_i` のビット幅)
- `0 <= sa, sb < 2^64`（LCG の初期 seed）

クエリは Knuth LCG（`s <- s * 6364136223846793005 + 1442695040888963407`）を 2 系列回し、
`a_i = sa_i & ((1 << abits) - 1)`, `b_i = sb_i & ((1 << abits) - 1)` として `N` 個生成する。

## 出力

`acc = XOR_{i=0..N-1} gcd(a_i, b_i)` を 1 行で出力する。

## 比較対象 algos

- `naive_mod`: 教科書的 `a %= b` Euclidean
- `std_gcd`: `std::gcd` (libstdc++ 実装)
- `binary`: Stein のアルゴリズム (binary GCD, 分岐あり)
- `binary_branchless`: Stein を `min/max` (cmov) で書いた変種
- `extgcd`: 拡張 Euclid (Bezout 係数も計算するわざと遅い参照点)
