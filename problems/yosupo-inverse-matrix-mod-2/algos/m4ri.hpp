#pragma once
#include "_common.hpp"
// M4RI で逆行列を計算。[A | I] (N×2N) を 8 列ブロックの Gauss-Jordan。
//
// 前進消去: 8 列ずつブロック内で k 個のピボット → TMP[256] でブロック外行を一括 XOR。
// 後退消去: 同じく 8 列ずつ、上の行を TMP で消去 (RREF にする)。
//
// 計算量: O(N³ / (64·k)) ≒ k=8 で plain の 1/8。
struct Inv {
 static vector<string> run(int n_in, const vector<string>& a) {
  constexpr int K = 8;
  const int N = n_in;
  const int W = (2 * N + 63) / 64;  // [A | I] の幅
  vector<u64> M((size_t) N * W, 0);
  for (int i = 0; i < N; ++i) {
   for (int j = 0; j < N; ++j)
    if (a[i][j] == '1') M[(size_t) i * W + j / 64] |= u64(1) << (j % 64);
   // 右半分の単位行列
   int col = N + i;
   M[(size_t) i * W + col / 64] |= u64(1) << (col % 64);
  }

  vector<u64> TMP((size_t) 256 * W, 0);

  // 前進消去 (cols [0, N) のみピボットとして使う; 右半分は副産物)
  for (int c = 0; c < N; c += K) {
   int k = std::min(K, N - c);
   for (int j = 0; j < k; ++j) {
    int target_col = c + j;
    int target_row = c + j;  // 前進では pivot_top = c (NxN 全部成功する前提)
    int sel = -1;
    for (int i = target_row; i < N; ++i) {
     if ((M[(size_t) i * W + target_col / 64] >> (target_col % 64)) & 1) { sel = i; break; }
    }
    if (sel == -1) return {};  // 可逆でない
    if (sel != target_row) {
     for (int w = 0; w < W; ++w) std::swap(M[(size_t) sel * W + w], M[(size_t) target_row * W + w]);
    }
    // batch 内の他行から target_col を消去
    for (int i = c; i < c + k; ++i) {
     if (i != target_row && ((M[(size_t) i * W + target_col / 64] >> (target_col % 64)) & 1)) {
      for (int w = 0; w < W; ++w) M[(size_t) i * W + w] ^= M[(size_t) target_row * W + w];
     }
    }
   }
   // TMP[2^k] 構築
   for (int w = 0; w < W; ++w) TMP[w] = 0;
   for (int s = 1; s < (1 << k); ++s) {
    int low = s & -s;
    int b = __builtin_ctz(low);
    int prev = s ^ low;
    for (int w = 0; w < W; ++w) TMP[(size_t) s * W + w] = TMP[(size_t) prev * W + w] ^ M[(size_t) (c + b) * W + w];
   }
   // batch 外の全行に適用 ([0, c) と [c+k, N))
   for (int i = 0; i < N; ++i) {
    if (i >= c && i < c + k) continue;
    int mask = (M[(size_t) i * W + c / 64] >> (c % 64)) & ((1 << k) - 1);
    if (mask) {
     for (int w = 0; w < W; ++w) M[(size_t) i * W + w] ^= TMP[(size_t) mask * W + w];
    }
   }
  }
  // この時点で左半分は単位行列、右半分は逆行列。
  vector<string> out(N, string(N, '0'));
  for (int i = 0; i < N; ++i) {
   for (int j = 0; j < N; ++j) {
    int col = N + j;
    if ((M[(size_t) i * W + col / 64] >> (col % 64)) & 1) out[i][j] = '1';
   }
  }
  return out;
 }
};
