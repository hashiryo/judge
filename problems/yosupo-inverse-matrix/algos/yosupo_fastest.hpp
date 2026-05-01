#pragma once
// =============================================================================
// Source: yosupo "Inverse Matrix" 提出 314161
//   https://judge.yosupo.jp/submission/314161
//   方式: Montgomery (R = 2^32) + AVX2 で flat[n*stride] と out[n*stride] (= I)
//   を持ち、Gauss-Jordan で out ← A^{-1} を構成。stride は 4 の倍数に
//   padding して _mm256_load_si256 / _mm256_store_si256 を直接使う。
//   1 行ぶんの mul-add (in Montgomery) を 4 lane で書き、減算は
//   `cmpgt` + `andnot` で実現。
// 抽出方針:
//   - I/O 部 (cin.read 高速 parser + cout.write) は base.cpp に移譲
//   - 行列値 (flat / out) は struct のローカル aligned_alloc で確保
//   - 元 main の return; (-1 出力) は InverseResult{false, {}} に置き換え
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx,avx2,sse")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#elif defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "_common.hpp"

namespace yosupo_314161 {
using ll = uint64_t;
constexpr ll P = 998244353;

constexpr ll multInv_const(ll x, ll m) {
 return 1 < (x %= m) ? m - multInv_const(m, x) * m / x : 1;
}
inline ll multInv(ll x, ll m) {
 return 1 < (x %= m) ? m - multInv(m, x) * m / x : 1;
}

inline InverseResult solve(int n, const vector<vector<u32>>& a_in) {
 constexpr ll R = 1ull << 32;
 constexpr ll R2 = (R % P) * R % P;
 constexpr ll R3 = R2 * R % P;
 const ll inv = multInv_const(R - P, R);
 auto divR = [](ll x) -> ll { return x >> 32; };
 auto modR = [](ll x) -> ll { return x & (R - 1); };
 auto redc = [&](ll v) -> ll {
  ll x = modR(modR(v) * inv);
  v = divR(v + x * P);
  if (v >= P) v -= P;
  return v;
 };
 auto mul = [&](ll a, ll b) -> ll { return redc(a * b); };

 int stride = (n + 3) & ~3;
 ll* flat = (ll*) std::aligned_alloc(32, (size_t) n * stride * 8);
 ll* out  = (ll*) std::aligned_alloc(32, (size_t) n * stride * 8);
 std::memset(flat, 0, (size_t) n * stride * 8);
 std::memset(out, 0, (size_t) n * stride * 8);
 for (int i = 0; i < n; ++i) {
  for (int j = 0; j < n; ++j) {
   flat[i * stride + j] = mul((ll) a_in[i][j], R2);
   out[i * stride + j]  = (i == j) ? (R % P) : 0;
  }
 }

 auto loop = [&](ll* src, ll* dst, int start, int end, ll f) -> void {
  __m256i invs = _mm256_set1_epi64x((long long) inv);
  __m256i mods = _mm256_set1_epi64x((long long) P);
  __m256i fs   = _mm256_set1_epi64x(f == 0 ? 0 : (long long) (P - f));
  for (int j = start & ~3; j < end; j += 4) {
   __m256i xr = _mm256_load_si256((__m256i*) (src + j));
   __m256i prod = _mm256_mul_epu32(fs, xr);
   __m256i x = _mm256_mul_epu32(prod, invs);
   x = _mm256_mul_epu32(x, mods);
   prod = _mm256_add_epi64(prod, x);
   prod = _mm256_srli_epi64(prod, 32);
   __m256i mask = _mm256_cmpgt_epi64(mods, prod);
   mask = _mm256_andnot_si256(mask, mods);
   prod = _mm256_sub_epi64(prod, mask);
   __m256i xi = _mm256_load_si256((__m256i*) (dst + j));
   prod = _mm256_add_epi64(prod, xi);
   mask = _mm256_cmpgt_epi64(mods, prod);
   mask = _mm256_andnot_si256(mask, mods);
   prod = _mm256_sub_epi64(prod, mask);
   _mm256_store_si256((__m256i*) (dst + j), prod);
  }
 };

 // Forward elimination
 for (int r = 0; r < n; ++r) {
  for (int i = r; i < n; ++i) {
   if (flat[i * stride + r] != 0) {
    std::ranges::swap_ranges(flat + r * stride, flat + r * stride + stride,
                             flat + i * stride, flat + i * stride + stride);
    std::ranges::swap_ranges(out + r * stride, out + r * stride + stride,
                             out + i * stride, out + i * stride + stride);
    break;
   }
  }
  if (flat[r * stride + r] == 0) {
   std::free(flat); std::free(out);
   return {false, {}};
  }
  ll f = mul(multInv(flat[r * stride + r], P), R3);
  for (int j = 0; j < stride; ++j) flat[r * stride + j] = mul(flat[r * stride + j], f);
  for (int j = 0; j < stride; ++j) out[r * stride + j]  = mul(out[r * stride + j], f);
  for (int i = r + 1; i < n; ++i) {
   f = flat[i * stride + r];
   loop(flat + r * stride, flat + i * stride, r, stride, f);
   loop(out  + r * stride, out  + i * stride, 0, r + 1, f);
  }
 }
 // Back substitution
 for (int r = n - 1; r >= 0; --r) {
  for (int i = 0; i < r; ++i) {
   ll f = flat[i * stride + r];
   flat[i * stride + r] = 0;
   loop(out + r * stride, out + i * stride, 0, stride, f);
  }
 }

 vector<vector<u32>> ans(n, vector<u32>(n));
 for (int i = 0; i < n; ++i)
  for (int j = 0; j < n; ++j)
   ans[i][j] = (u32) redc(out[i * stride + j]);
 std::free(flat); std::free(out);
 return {true, std::move(ans)};
}
} // namespace yosupo_314161

struct Inverse {
 static InverseResult run(int n, const vector<vector<u32>>& a) {
  return yosupo_314161::solve(n, a);
 }
};
