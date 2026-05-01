#pragma once
#include "_common.hpp"
// M4RI でランクを計算。N < M なら転置 (rank 不変、cols 小さい方が速い)。
//
// rank では列にピボットが無い場合があるので、8 列固定ブロックではなく
// 「次の K 個のピボット列を batch として収集」する形で M4RI を適用。
// batch 内の列は不連続なことがあるので mask 抽出は各 col bit を 1 つずつ拾う。
struct Rank {
 static int run(int n_in, int m_in, const vector<string>& a) {
  constexpr int K = 8;
  bool transposed = (n_in < m_in);
  int rows = transposed ? m_in : n_in;
  int cols = transposed ? n_in : m_in;
  const int W = (cols + 63) / 64;
  vector<u64> M((size_t) rows * W, 0);
  if (transposed) {
   for (int i = 0; i < n_in; ++i)
    for (int j = 0; j < m_in; ++j)
     if (a[i][j] == '1') M[(size_t) j * W + i / 64] |= u64(1) << (i % 64);
  } else {
   for (int i = 0; i < n_in; ++i)
    for (int j = 0; j < m_in; ++j)
     if (a[i][j] == '1') M[(size_t) i * W + j / 64] |= u64(1) << (j % 64);
  }

  vector<u64> TMP((size_t) 256 * W, 0);
  int pivot_top = 0;
  int c = 0;
  while (c < cols && pivot_top < rows) {
   int batch_cols[K];
   int batch_size = 0;
   int cc = c;
   // K_eff = この batch で扱える最大行数 (実 row 数で頭打ち)
   const int K_eff = std::min(K, rows - pivot_top);
   while (batch_size < K_eff && cc < cols) {
    int target_row = pivot_top + batch_size;
    int sel = -1;
    for (int i = target_row; i < rows; ++i) {
     if ((M[(size_t) i * W + cc / 64] >> (cc % 64)) & 1) { sel = i; break; }
    }
    if (sel == -1) { cc++; continue; }
    if (sel != target_row) {
     for (int w = 0; w < W; ++w) std::swap(M[(size_t) sel * W + w], M[(size_t) target_row * W + w]);
    }
    // 重要: in-batch elim は pt 以下の全行 [pt, rows) を対象にする。
    // 範囲を [pt, pt+K_eff) に絞ると、後の batch col で「外から swap してきた
    // pivot 行が earlier batch col のビットを汚染した状態」で他行に XOR して
    // しまい、既存のピボット行の pivot ビットを破壊するバグになる。
    // 全行 elim にすれば外から swap してきた行は事前に earlier batch col 全部
    // クリーン済みなので安全。M4RI の "k 倍速化" は失われるが、これは
    // upper-triangular Gauss では本質的に避けられない (M4RI は本来 RREF
    // 後段の clean up を加速するもの)。
    for (int i = pivot_top; i < rows; ++i) {
     if (i != target_row && ((M[(size_t) i * W + cc / 64] >> (cc % 64)) & 1)) {
      for (int w = 0; w < W; ++w) M[(size_t) i * W + w] ^= M[(size_t) target_row * W + w];
     }
    }
    batch_cols[batch_size] = cc;
    batch_size++;
    cc++;
   }
   if (batch_size == 0) break;
   // TMP[2^batch_size] = batch ピボット行の XOR 部分集合
   for (int w = 0; w < W; ++w) TMP[w] = 0;
   for (int s = 1; s < (1 << batch_size); ++s) {
    int low = s & -s;
    int b = __builtin_ctz(low);
    int prev = s ^ low;
    for (int w = 0; w < W; ++w) TMP[(size_t) s * W + w] = TMP[(size_t) prev * W + w] ^ M[(size_t) (pivot_top + b) * W + w];
   }
   // batch_cols が連続なら mask 抽出を 1 byte read で
   bool consecutive = true;
   for (int j = 1; j < batch_size; ++j)
    if (batch_cols[j] != batch_cols[0] + j) { consecutive = false; break; }
   for (int i = 0; i < rows; ++i) {
    if (i >= pivot_top && i < pivot_top + batch_size) continue;
    int mask;
    if (consecutive) {
     mask = (int) ((M[(size_t) i * W + batch_cols[0] / 64] >> (batch_cols[0] % 64)) & ((1 << batch_size) - 1));
    } else {
     mask = 0;
     for (int j = 0; j < batch_size; ++j) {
      int col = batch_cols[j];
      mask |= (int) ((M[(size_t) i * W + col / 64] >> (col % 64)) & 1) << j;
     }
    }
    if (mask) {
     for (int w = 0; w < W; ++w) M[(size_t) i * W + w] ^= TMP[(size_t) mask * W + w];
    }
   }
   pivot_top += batch_size;
   c = cc;
  }
  return pivot_top;
 }
};
