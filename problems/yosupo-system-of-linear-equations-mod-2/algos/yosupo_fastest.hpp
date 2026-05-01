#pragma once
// =============================================================================
// Source: yosupo "System of Linear Equations (mod 2)" 提出 268155 を抽出。
//   https://judge.yosupo.jp/submission/268155
//   方式:
//   1. 拡張行列 [A | b] を bitset で持って RREF (前進消去のみ) する
//   2. 不整合 (0...0|1 行) があれば -1
//   3. ピボット列ごとに rref を再配置、自由変数 (FV[i] != 0) を抽出
//   4. 後ろから 8 列ずつブロック処理、Four Russians (256 entries TMP table) で
//      pivot 行への反映を加速 (M4RI 風の back-substitution)
//   5. sols[i] (M 個の bitset) の bit j を取って (R+1) × M の解行列を出力
// 抽出方針:
//   - I/O は base.cpp 経由
//   - bitset<4097> 周りの reinterpret_cast はそのまま
//   - 全体を namespace yosupo_268155 で囲む
//   - Solve::run wrapper で string ↔ bitset 変換 + 出力整形
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#pragma GCC optimize("Ofast,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#elif defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "_common.hpp"

namespace yosupo_268155 {
constexpr int SZ = 4097;
using AsArray = std::array<long, 65>;

inline std::vector<std::bitset<SZ>> RREF(std::vector<std::bitset<SZ>> x, int N, int M) {
 int h = 0, k = 0;
 while (h < N && k < M) {
  int idx = h;
  while (idx < (int) x.size() && x[idx][k] == 0) idx++;
  if (idx != (int) x.size()) {
   std::swap(x[h], x[idx]);
   for (int i = h + 1; i < (int) x.size(); i++)
    if (x[i][k]) x[i] ^= x[h];
   h++;
  }
  k++;
 }
 return x;
}
} // namespace yosupo_268155

struct Solve {
 static vector<string> run(int N, int M, const vector<string>& a, const vector<int>& bvec) {
  using namespace yosupo_268155;
  std::vector<std::bitset<SZ>> raw(std::max(N, M));
  for (int i = 0; i < N; ++i) {
   for (int j = 0; j < M; ++j)
    if (a[i][j] == '1') raw[i][j] = 1;
   raw[i][M] = bvec[i] & 1;
  }
  auto rref = RREF(raw, N, M);
  // 不整合チェック: count() == 1 かつ M 列目だけ立つ → 解なし
  for (int i = 0; i < std::max(N, M); ++i)
   if (rref[i].count() == 1 && rref[i][M]) return {};
  // 非零行を集める
  std::vector<std::bitset<SZ>> pruned_rref;
  for (int i = 0; i < std::max(N, M); ++i)
   if (rref[i].count()) pruned_rref.push_back(rref[i]);
  // 各行の lead 列を見つけて rref[lead] にセット (lead で sort 済の状態に再配置)
  rref = std::vector<std::bitset<SZ>>(M);
  for (int i = 0; i < (int) pruned_rref.size(); ++i) {
   int lead = 0;
   while (lead < M && !pruned_rref[i][lead]) lead++;
   rref[lead] = pruned_rref[i];
  }
  // 自由変数の特定: rref[i][i] == 0 なら i は自由変数。FV[i] = (1-based index)。
  std::vector<int> fv, FV(M);
  for (int i = 0; i < M; ++i)
   if (!rref[i][i]) {
    FV[i] = (int) fv.size() + 1;
    fv.push_back(i);
   }
  int R = (int) fv.size();
  int MM = (M + 7) / 8 * 8;
  std::vector<std::bitset<SZ>> sols(MM);
  std::vector<std::bitset<SZ>> TMP(256);
  for (int first = MM - 8; first >= 0; first -= 8) {
   for (int i = std::min(M - 1, first + 7); i >= first; --i) {
    if (FV[i]) sols[i][FV[i]] = 1;
    else {
     sols[i][0] = sols[i][0] ^ rref[i][M];
     for (int j = i + 1; j < std::min(M, first + 8); ++j)
      if (rref[i][j]) sols[i] ^= sols[j];
    }
   }
   // M4RI: 8-bit 部分集合 XOR テーブル
   for (int bit = 0; bit < 8; ++bit)
    for (int cur = 0; cur < (1 << bit); ++cur)
     TMP[cur + (1 << bit)] = TMP[cur] ^ sols[first + bit];
   for (int i = 0; i < std::min(first, M); ++i) {
    if (!FV[i]) {
     AsArray& aa = *reinterpret_cast<AsArray*>(&rref[i]);
     const int mask = (aa[first >> 6] >> (first & 63)) & 255;
     sols[i] ^= TMP[mask];
    }
   }
  }
  // sols[i] の bit j (j=0..R) を (R+1) × M の解行列に展開: out[j][i] = sols[i][j]
  vector<string> out(R + 1, string(M, '0'));
  for (int i = 0; i < M; ++i) {
   AsArray& aa = *reinterpret_cast<AsArray*>(&sols[i]);
   for (int j = 0; j <= R; ++j)
    if ((aa[j >> 6] >> (j & 63)) & 1) out[j][i] = '1';
  }
  return out;
 }
};
