#pragma once
// =============================================================================
// Source: yosupo "Matrix Product (mod 998244353)" 提出 193471 を抽出。
//   https://judge.yosupo.jp/submission/193471
//   方式: AVX2 + Montgomery reduction で 64×64 base kernel を組み、
//   その上に再帰 Strassen (7-mul) を載せた行列積。Morton 順に再配置して
//   キャッシュ効率を稼ぐ + Strassen で N=1024 で N³ → N^log2(7) ≒ N^2.81 削減。
// 抽出方針:
//   - I/O 部分 (Qinf/Qoutf) は base.cpp に移譲して削除
//   - 全体を namespace yosupo_193471 に閉じる
//   - 静的 global a/b/c/A/B/C 配列はそのまま namespace 内
//   - Mul::run wrapper で flat vector<u32> ↔ 静的 array (パディング込み)
// 制約: N, M, P ≤ 1024 を前提。
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#elif defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "_common.hpp"

namespace yosupo_193471 {
using i64 = int64_t;
using u8 = uint8_t;
using u16 = uint16_t;
using u32_ = uint32_t;
using u64_ = uint64_t;
using I256 = __m256i;
using idt = std::size_t;

constexpr u32_ mod = 998244353;
constexpr u32_ niv = [] { u32_ n = 2 + mod; for (int i = 0; i < 4; ++i) { n*= 2 + mod * n; } return n; }();
constexpr u32_ R2 = (u32_) (-(u64_) mod % mod);
constexpr i64 ngm = i64(1) << 63;

constexpr u32_ reduce_s(u64_ x, u32_ niv_, u32_ M) { return (x + (u64_) (u32_(x) * niv_) * M) >> 32; }
constexpr u32_ shrk32_s(u32_ x, u32_ M) { return std::min(x, x - M); }
inline I256 reduce_v(I256 a, I256 b, I256 niv_, I256 M) {
 I256 kil = _mm256_mul_epu32(a, niv_), jok = _mm256_mul_epu32(b, niv_);
 kil = _mm256_mul_epu32(kil, M), jok = _mm256_mul_epu32(jok, M);
 return _mm256_blend_epi32(_mm256_srli_epi64(_mm256_add_epi64(a, kil), 32), _mm256_add_epi64(b, jok), 0xaa);
}
inline I256 mul_sm(I256 a, I256 b, I256 niv_, I256 M) {
 return reduce_v(_mm256_mul_epu32(a, b), _mm256_mul_epu32(_mm256_srli_epi64(a, 32), b), niv_, M);
}
inline I256 shrk32(I256 x, I256 M) { return _mm256_min_epu32(x, _mm256_sub_epi32(x, M)); }
inline I256 dilt32(I256 x, I256 M) { return _mm256_min_epu32(x, _mm256_add_epi32(x, M)); }

constexpr idt BLK = 64;

[[gnu::noinline]] void __kernel_1(const u32_* __restrict__ a, const u32_* __restrict__ b, u32_* __restrict__ c) {
 const I256 M = _mm256_set1_epi32((int) mod), M2 = _mm256_set1_epi32((int) (mod * 2));
 const I256 Niv = _mm256_set1_epi32((int) niv), Ngm = _mm256_set1_epi64x(ngm);
#define INIT_R(p) I256 r##p = _mm256_setzero_si256(), R##p = _mm256_setzero_si256();
#define YZLF_WORK(p) r##p = _mm256_add_epi64(r##p, _mm256_mul_epu32(_mm256_set1_epi32((int) a[(i + p) * BLK + w]), z0)), R##p = _mm256_add_epi64(R##p, _mm256_mul_epu32(_mm256_set1_epi32((int) a[(i + p) * BLK + w]), Z0));
#define R_SHRK(p) r##p = _mm256_sub_epi64(r##p, _mm256_min_epu32(M2, _mm256_and_si256(Ngm, r##p))), R##p = _mm256_sub_epi64(R##p, _mm256_min_epu32(M2, _mm256_and_si256(Ngm, R##p)));
#define REDUCE_OP(p) _mm256_store_si256((I256*) (c + (i + p) * BLK + j), shrk32(shrk32(reduce_v(r##p, R##p, Niv, M), M2), M));
#define UNR8(d) d(0) d(1) d(2) d(3) d(4) d(5) d(6) d(7)
 for (idt i = 0; i < BLK; i+= 8) {
  for (idt j = 0; j < BLK; j+= 8) {
   UNR8(INIT_R)
   for (idt k = 0; k < BLK; k+= 8) {
    for (idt w = k; w < k + 8; ++w) {
     const I256 z0 = _mm256_load_si256((I256*) (b + w * BLK + j)), Z0 = _mm256_srli_epi64(z0, 32);
     UNR8(YZLF_WORK)
    }
    UNR8(R_SHRK)
   }
   UNR8(REDUCE_OP)
  }
 }
#undef INIT_R
#undef YZLF_WORK
#undef R_SHRK
#undef REDUCE_OP
#undef UNR8
}

// 静的 workspace。N=1024 で n²=1M, 4n²/3≒1.4M。
constexpr idt PAD_NN_MAX = 1024;
alignas(32) inline unsigned int A[(1 << 22) / 3];
alignas(32) inline unsigned int B[(1 << 22) / 3];
alignas(32) inline unsigned int C_buf[(1 << 22) / 3];
alignas(32) inline unsigned int a_buf[1 << 20];
alignas(32) inline unsigned int b_buf[1 << 20];
alignas(32) inline unsigned int c_buf[1 << 20];

// libc++ には std::__lg が無いので __builtin_clzll で代替。
constexpr idt bcl(idt x) { return x < 2 ? 1 : idt(2) << (63 - __builtin_clzll((unsigned long long) (x - 1))); }

inline void __place_mat(idt x, idt y, idt n, idt N, const u32_* __restrict__ a, u32_* __restrict__ Aout) {
 if (n == BLK) {
  for (idt i = 0; i < BLK; ++i) std::memcpy(Aout + i * BLK, a + (y + i) * N + x, BLK * sizeof(u32_));
  return;
 }
 idt nn = n / 2, D = nn * nn;
 __place_mat(x, y, nn, N, a, Aout);
 __place_mat(x + nn, y, nn, N, a, Aout + D);
 __place_mat(x, y + nn, nn, N, a, Aout + D * 2);
 __place_mat(x + nn, y + nn, nn, N, a, Aout + D * 3);
}

inline void __place_mat_anti_fx(idt x, idt y, idt n, idt N, u32_* __restrict__ a, const u32_* __restrict__ Ain) {
 if (n == BLK) {
  const I256 R2x8 = _mm256_set1_epi32((int) R2), Niv = _mm256_set1_epi32((int) niv), M = _mm256_set1_epi32((int) mod);
  for (idt i = 0; i < BLK; ++i) {
   for (idt j = 0; j < BLK; j+= 8) {
    _mm256_store_si256((I256*) (a + (y + i) * N + x + j),
     shrk32(mul_sm(_mm256_load_si256((const I256*) (Ain + i * BLK + j)), R2x8, Niv, M), M));
   }
  }
  return;
 }
 idt nn = n / 2, D = nn * nn;
 __place_mat_anti_fx(x, y, nn, N, a, Ain);
 __place_mat_anti_fx(x + nn, y, nn, N, a, Ain + D);
 __place_mat_anti_fx(x, y + nn, nn, N, a, Ain + D * 2);
 __place_mat_anti_fx(x + nn, y + nn, nn, N, a, Ain + D * 3);
}

inline void __Strassen_matmul(u32_* __restrict__ a, u32_* __restrict__ b, u32_* __restrict__ c, idt n,
                              u32_* __restrict__ aa, u32_* __restrict__ bb) {
 if (n == BLK) { __kernel_1(a, b, c); return; }
 idt nn = n / 2, D = nn * nn;
 u32_* cc = c + n * n;
 idt i00 = 0, i01 = D, i10 = D * 2, i11 = D * 3;
 auto add = [&](u32_* x, const u32_* y, const u32_* z) {
  const I256 M = _mm256_set1_epi32((int) mod);
  for (idt i = 0; i < D; i+= 8) {
   _mm256_store_si256((I256*) (x + i),
    shrk32(_mm256_add_epi32(_mm256_load_si256((const I256*) (y + i)), _mm256_load_si256((const I256*) (z + i))), M));
  }
 };
 auto sub = [&](u32_* x, const u32_* y, const u32_* z) {
  const I256 M = _mm256_set1_epi32((int) mod);
  for (idt i = 0; i < D; i+= 8) {
   _mm256_store_si256((I256*) (x + i),
    dilt32(_mm256_sub_epi32(_mm256_load_si256((const I256*) (y + i)), _mm256_load_si256((const I256*) (z + i))), M));
  }
 };
 sub(aa, a + i01, a + i11);
 add(bb, b + i10, b + i11);
 __Strassen_matmul(aa, bb, c + i00, nn, aa + D, bb + D);

 add(aa, a + i00, a + i01);
 __Strassen_matmul(aa, b + i11, c + i01, nn, aa + D, bb + D);
 sub(c + i00, c + i00, c + i01);

 sub(bb, b + i10, b + i00);
 __Strassen_matmul(a + i11, bb, c + i10, nn, aa + D, bb + D);
 add(c + i00, c + i00, c + i10);

 sub(aa, a + i10, a + i00);
 add(bb, b + i00, b + i01);
 __Strassen_matmul(aa, bb, c + i11, nn, aa + D, bb + D);

 add(aa, a + i10, a + i11);
 __Strassen_matmul(aa, b + i00, cc, nn, aa + D, bb + D);
 add(c + i10, c + i10, cc);
 sub(c + i11, c + i11, cc);

 add(aa, a + i00, a + i11);
 add(bb, b + i00, b + i11);
 __Strassen_matmul(aa, bb, cc, nn, aa + D, bb + D);
 add(c + i00, c + i00, cc);
 add(c + i11, c + i11, cc);

 sub(bb, b + i01, b + i11);
 __Strassen_matmul(a + i00, bb, cc, nn, aa + D, bb + D);
 add(c + i01, c + i01, cc);
 add(c + i11, c + i11, cc);
}

inline void mul_Strassen(const u32_* __restrict__ a, const u32_* __restrict__ b, u32_* __restrict__ c, idt N) {
 __place_mat(0, 0, N, N, a, A);
 __place_mat(0, 0, N, N, b, B);
 __Strassen_matmul(A, B, C_buf, N, A + N * N, B + N * N);
 __place_mat_anti_fx(0, 0, N, N, c, C_buf);
}

} // namespace yosupo_193471

struct Mul {
 static vector<u32> run(int n, int m, int p, const vector<u32>& a, const vector<u32>& b) {
  using namespace yosupo_193471;
  idt N = std::max(BLK, bcl((idt) std::max({n, m, p})));
  // padded square N×N に詰める。未使用領域は 0。
  std::memset(a_buf, 0, sizeof(a_buf));
  std::memset(b_buf, 0, sizeof(b_buf));
  std::memset(c_buf, 0, sizeof(c_buf));
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < m; ++j) a_buf[i * N + j] = a[(size_t) i * m + j];
  for (int i = 0; i < m; ++i)
   for (int j = 0; j < p; ++j) b_buf[i * N + j] = b[(size_t) i * p + j];
  mul_Strassen(a_buf, b_buf, c_buf, N);
  vector<u32> result((size_t) n * p);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < p; ++j) result[(size_t) i * p + j] = c_buf[i * N + j];
  return result;
 }
};
