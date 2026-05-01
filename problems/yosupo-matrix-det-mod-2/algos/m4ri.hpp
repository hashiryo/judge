#pragma once
#include "_common.hpp"
// Method of Four Russians (M4RI) で det を計算。
// 8 列ずつブロック処理 + 256 entries の XOR テーブル参照で消去を加速。
//
// アイデア:
//   1. 8 列の中で k 個 (= 8) のピボットを見つけ、そのピボット行内だけ
//      互いの 8 列を消去 (k×k の部分 Gauss)。pivot 行数: O(k²·W)
//   2. TMP[mask] = mask の各 bit が 1 のピボット行を XOR したもの。256 entries。
//   3. 他の全ての行 i について、cols [c, c+8) の 8 bit を mask として
//      row[i] ^= TMP[mask] を 1 回 XOR。
//
// 計算量: O(N³ / (64·k)) で plain Gauss の O(N³/64) より定数倍 (≒8x) 速い。
//
// det 専用なので、ピボットが見つからなければ即 return 0 (行列特異 = det 0)。
struct Det {
 static int run(int n, const vector<string>& a) {
  constexpr int K = 8;
  const int W = (n + 63) / 64;
  vector<u64> M((size_t) n * W, 0);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j)
    if (a[i][j] == '1') M[(size_t) i * W + j / 64] |= u64(1) << (j % 64);

  vector<u64> TMP((size_t) 256 * W, 0);
  int pivot_top = 0;
  for (int c = 0; c < n; c += K) {
   int k = std::min(K, n - c);
   // Step 1: ブロック内 Gauss で k 個のピボットを見つける
   for (int j = 0; j < k; ++j) {
    int target_col = c + j;
    int target_row = pivot_top + j;
    int sel = -1;
    for (int i = target_row; i < n; ++i) {
     if ((M[(size_t) i * W + target_col / 64] >> (target_col % 64)) & 1) { sel = i; break; }
    }
    if (sel == -1) return 0;  // 特異行列 → det = 0
    if (sel != target_row) {
     for (int w = 0; w < W; ++w) std::swap(M[(size_t) sel * W + w], M[(size_t) target_row * W + w]);
    }
    // 重要: in-batch elim は全行 [pt, n) を対象にする。
    // 範囲を [pt, pt+k) に絞ると、後の batch col で外部から swap してきた
    // pivot 行が earlier batch col のビットを汚染した状態で他行に XOR してしまい、
    // 既存のピボット行の pivot ビットを破壊するバグになる。
    for (int i = pivot_top; i < n; ++i) {
     if (i != target_row && ((M[(size_t) i * W + target_col / 64] >> (target_col % 64)) & 1)) {
      for (int w = 0; w < W; ++w) M[(size_t) i * W + w] ^= M[(size_t) target_row * W + w];
     }
    }
   }
   // Step 2: TMP[2^k] を構築。TMP[s] = bit 立てた pivot 行の XOR。
   for (int w = 0; w < W; ++w) TMP[w] = 0;
   for (int s = 1; s < (1 << k); ++s) {
    int low = s & -s;
    int b = __builtin_ctz(low);
    int prev = s ^ low;
    for (int w = 0; w < W; ++w) TMP[(size_t) s * W + w] = TMP[(size_t) prev * W + w] ^ M[(size_t) (pivot_top + b) * W + w];
   }
   // Step 3: pivot 行以外の全行に TMP を適用
   for (int i = 0; i < n; ++i) {
    if (i >= pivot_top && i < pivot_top + k) continue;
    int mask = (M[(size_t) i * W + c / 64] >> (c % 64)) & ((1 << k) - 1);
    if (mask) {
     for (int w = 0; w < W; ++w) M[(size_t) i * W + w] ^= TMP[(size_t) mask * W + w];
    }
   }
   pivot_top += k;
  }
  return 1;
 }
};
