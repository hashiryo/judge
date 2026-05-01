#pragma once
#include "_common.hpp"
// Four Russians (M4RI 風): 8 行ずつ「全部分集合の XOR テーブル」を作って参照。
// O(N*M*K / log) 級。bitset の素朴版より定数倍速い。
struct Mul {
 static constexpr int MAXN = 4096;
 static vector<string> run(int n, int m, int k, const vector<string>& a, const vector<string>& b) {
  // B を u64 行 (各行 K bits パック) で持つ。
  int kw = (k + 63) / 64;
  vector<u64> bw((size_t) m * kw, 0);
  for (int i = 0; i < m; ++i)
   for (int j = 0; j < k; ++j)
    if (b[i][j] == '1') bw[(size_t) i * kw + j / 64] |= u64(1) << (j % 64);
  // A も同様に u64 行 (M bits パック) として持つが、Four Russians のキー抽出のために
  // 各 8 ビット (バイト) ごとに 256 エントリのテーブルを引いて結果を XOR していく。
  int mw8 = (m + 7) / 8;
  vector<u8> ab((size_t) n * mw8, 0);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < m; ++j)
    if (a[i][j] == '1') ab[(size_t) i * mw8 + j / 8] |= u8(1) << (j % 8);
  // テーブル: 各 8 行ブロックについて、ビットの全 256 部分集合の XOR を bw 行から作る。
  vector<u64> table((size_t) 256 * kw, 0);
  vector<u64> res((size_t) n * kw, 0);
  for (int blk = 0; blk < m; blk += 8) {
   int rows = std::min(8, m - blk);
   // table[mask] = XOR of bw[blk + bit] for each bit set in mask
   std::fill(table.begin(), table.end(), 0);
   for (int s = 1; s < (1 << rows); ++s) {
    int low = s & -s;
    int p = __builtin_ctz(low);
    int prev = s ^ low;
    for (int j = 0; j < kw; ++j) {
     table[(size_t) s * kw + j] = table[(size_t) prev * kw + j] ^ bw[(size_t) (blk + p) * kw + j];
    }
   }
   for (int i = 0; i < n; ++i) {
    int mask = ab[(size_t) i * mw8 + blk / 8];
    if (rows < 8) mask &= (1 << rows) - 1;
    if (!mask) continue;
    for (int j = 0; j < kw; ++j) res[(size_t) i * kw + j] ^= table[(size_t) mask * kw + j];
   }
  }
  vector<string> out(n, string(k, '0'));
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < k; ++j)
    if (res[(size_t) i * kw + j / 64] >> (j % 64) & 1) out[i][j] = '1';
  return out;
 }
};
