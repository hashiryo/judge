# warshall-floyd

全点対最短路 (Warshall-Floyd) の SIMD 比較ベンチ問題。

距離行列を入力 seed から LCG で生成し、WF を回したあとの全 dist[i][j] (i, j ∈ [0, V)) の XOR を出す。

## 入力

`V W seed`

- `1 <= V`（頂点数）
- `1 <= W`（辺重みの上限。各辺重みは `[1, W]` 一様）
- `0 <= seed < 2^64`

dist は i32 で持っているので、`V * W < WF_INF = 10^9` を満たす範囲で振る
(テストケースは max_path ≤ 10^8 程度に収まるよう構成)。

`d[i][j]` 初期値は LCG 生成の `[1, W]`。`d[i][i] = 0`。

## 出力

WF 後の全 `d[i][j]` の XOR (u32 として、u64 にゼロ拡張して累積) を 1 行で出力する。

## 比較対象 algos

- `naive`: scalar 3 重ループ
- `simd`: AVX2 で内側 j ループを 8-wide ベクトル化 (`_mm256_min_epi32` + `_mm256_add_epi32`)

`Vp = ceil(V / 8) * 8` でパディングし、SIMD 版は j ループを Vp までベクトル化する
(パディング列は INF で埋まっており、min で影響しない)。
