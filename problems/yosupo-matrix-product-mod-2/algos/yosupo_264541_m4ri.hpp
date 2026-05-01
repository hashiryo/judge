#pragma once
#include "_common.hpp"
// =============================================================================
// Source: yosupo "Matrix Product (mod 2)" 提出 264541 を抽出。
//   https://judge.yosupo.jp/submission/264541
//   方式: Four Russians (M4RI) + 256-行ブロック + 8-bit mask テーブル参照。
//   提出は I/O も SIMD 高速化していたが、ここでは algo 部分のみ抽出 (I/O は base.cpp)。
//   N, M, K ≤ 4096 を前提とし、4096×4096 で固定処理。未使用領域は 0 で埋まる。
// 抽出方針:
//   - I/O は base.cpp 経由 (string ↔ u64 packing は内部で実施)
//   - 静的 global 配列 a/b/c/TMP は namespace 内に閉じる
//   - Mul::run() wrapper で N×M / M×K → N×K の string 変換
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2")
#endif

namespace yosupo_264541 {
constexpr int FN = 4096;
alignas(32) inline u64 a[FN][64];
alignas(32) inline u64 b[FN][64];
alignas(32) inline u64 c[FN][64];
alignas(32) inline u64 TMP[256][64];
} // namespace yosupo_264541

struct Mul {
 static vector<string> run(int n, int m, int k, const vector<string>& av, const vector<string>& bv) {
  using namespace yosupo_264541;
  std::memset(a, 0, sizeof(a));
  std::memset(b, 0, sizeof(b));
  std::memset(c, 0, sizeof(c));
  // a[i][j>>6] の (j&63) ビット目に av[i][j] を詰める。
  for (int i = 0; i < n; ++i) {
   const char* row = av[i].data();
   for (int j = 0; j < m; ++j)
    if (row[j] == '1') a[i][j >> 6] |= u64(1) << (j & 63);
  }
  for (int i = 0; i < m; ++i) {
   const char* row = bv[i].data();
   for (int j = 0; j < k; ++j)
    if (row[j] == '1') b[i][j >> 6] |= u64(1) << (j & 63);
  }
  // Four Russians: B の 8 行ブロックごとに TMP テーブル (256 entries) を構築、
  // A の各 8-bit mask に対応する TMP を XOR で C に積む。
  for (int j = 0; j < FN; j += 8) {
   // TMP[0] は 0 のまま。TMP[s + (1<<bit)] = TMP[s] ^ b[j+bit] を bit ごとに展開。
   std::memset(TMP[0], 0, sizeof(TMP[0]));
   for (int bit = 0; bit < 8; ++bit)
    for (int cur = 0; cur < (1 << bit); ++cur)
     for (int kw = 0; kw < 64; ++kw) TMP[cur + (1 << bit)][kw] = TMP[cur][kw] ^ b[j + bit][kw];
   for (int i = 0; i < FN; i += 256) {
    for (int X = 0; X < 256; ++X) {
     int mask = (a[i + X][j >> 6] >> (j & 63)) & 255;
     if (mask) {
      for (int kw = 0; kw < 64; ++kw) c[i + X][kw] ^= TMP[mask][kw];
     }
    }
   }
  }
  vector<string> res(n, string(k, '0'));
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < k; ++j)
    if ((c[i][j >> 6] >> (j & 63)) & 1) res[i][j] = '1';
  }
  return res;
 }
};
