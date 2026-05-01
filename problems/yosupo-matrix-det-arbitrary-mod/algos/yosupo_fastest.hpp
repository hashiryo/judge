#pragma once
// =============================================================================
// Source: yosupo "Matrix Determinant (arbitrary mod)" 提出 193004 を抽出 (1 位)。
//   https://judge.yosupo.jp/submission/193004
//   方式: AVX-512F + AVX-512DQ で 16 lane × u32 を一気に Barrett-like 乗算。
//   `fmd_32` (= floating-point Montgomery-like reduction): 1/M を double で
//   保持 (`nextafter` で safe round)、`u64*u64 → double で割って商を求めて
//   引き戻す」スタイル。AVX-512 で u64↔double 変換が直接できるのが効く。
//
//   x86_64 native: `#pragma GCC target("avx512f,avx512dq")` で AVX-512 命令
//   有効化 (CPU が対応していなければ実行時 SIGILL → RE で記録)。
//   ARM (USE_SIMDE): simde の AVX-512 emulation 経由で動くが、NEON は 128-bit
//   なので 4 命令/AVX-512 になり激遅 (動作するかも怪しい)。
// 抽出方針:
//   - I/O は base.cpp 経由
//   - 全体を namespace yosupo_193004 で囲む
//   - Det::run wrapper で N, mod, A → u32 結果
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
// AVX-512F + AVX-512DQ が必要。simde の AVX-512DQ サポートは
// 不完全 (_mm512_cvt_roundepu64_pd 等が未実装) なので x86_64 native 限定。
// 非 x86_64 環境では runtime abort の stub を提供。
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx512f,avx512dq")
#include <immintrin.h>
#define YOSUPO_193004_ENABLE 1
#endif

#include "_common.hpp"

#ifndef YOSUPO_193004_ENABLE
// 非 x86_64 native 環境: runtime abort
struct Det {
 static u32 run(int, u32, const vector<vector<u32>>&) {
  fprintf(stderr, "yosupo_193004 (AVX-512): x86_64 + AVX-512 required\n");
  std::abort();
 }
};
#else

namespace yosupo_193004 {
using f64 = double;
using vd8 = __m512d;
using I512 = __m512i;

struct fmd_32 {
 f64 rp;
 uint32_t M;
 fmd_32(uint32_t m): rp(std::nextafter(1.0 / m, 1.0)), M(m) {}
 uint32_t operator()(uint64_t x) const {
  uint32_t res = (uint32_t) (x - (uint64_t) (f64(uint64_t(x)) * rp) * M);
  return std::min(res, res + M);
 }
};

struct fmd_32_I512 {
 vd8 rp;
 I512 M;
 fmd_32_I512(fmd_32 m): rp(_mm512_set1_pd(m.rp)), M(_mm512_set1_epi32((int) m.M)) {}
 // x*y+z (Barrett-like via floating-point)
 I512 mul(I512 x, I512 y, I512 z) const {
  I512 a = _mm512_mul_epu32(x, y);
  I512 b = _mm512_mul_epu32(_mm512_srli_epi64(x, 32), _mm512_srli_epi64(y, 32));
  a = _mm512_add_epi64(a, _mm512_srli_epi64(_mm512_slli_epi64(z, 32), 32));
  b = _mm512_add_epi64(b, _mm512_srli_epi64(z, 32));
  vd8 c = _mm512_cvt_roundepu64_pd(a, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  vd8 d = _mm512_cvt_roundepu64_pd(b, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  c = _mm512_mul_pd(c, rp);
  d = _mm512_mul_pd(d, rp);
  I512 e = _mm512_cvt_roundpd_epu64(c, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  I512 f = _mm512_cvt_roundpd_epu64(d, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
  a = _mm512_sub_epi64(a, _mm512_mul_epu32(e, M));
  b = _mm512_sub_epi64(b, _mm512_mul_epu32(f, M));
  I512 res = _mm512_mask_blend_epi32(0xaaaa, a, _mm512_slli_epi64(b, 32));
  return _mm512_min_epu32(res, _mm512_add_epi32(res, M));
 }
};

inline uint32_t det_impl(vector<vector<uint32_t>> A, const fmd_32 fm) {
 uint32_t res = 1;
 int n = (int) A.size();
 const fmd_32_I512 fms(fm);
 for (int i = 0; i < n; ++i) {
  for (int j = i + 1; j < n; std::swap(A[i], A[j]), res = fm.M - res, ++j) {
   for (; A[i][i]; std::swap(A[i], A[j]), res = fm.M - res) {
    uint32_t o = fm.M - (A[j][i] / A[i][i]);
    int k = i;
    if (n - k >= 16) {
     I512 ox16 = _mm512_set1_epi32((int) o);
     for (; k + 15 < n; k += 16) {
      _mm512_storeu_si512((I512*) (A[j].data() + k),
       fms.mul(_mm512_loadu_si512((const I512*) (A[i].data() + k)), ox16,
        _mm512_loadu_si512((const I512*) (A[j].data() + k))));
     }
    }
    for (; k < n; ++k) {
     A[j][k] = fm(A[j][k] + uint64_t(o) * A[i][k]);
    }
   }
  }
  res = fm(res * uint64_t(A[i][i]));
  if (!res) break;
 }
 return res;
}
} // namespace yosupo_193004

struct Det {
 static u32 run(int n, u32 mod, const vector<vector<u32>>& a) {
  yosupo_193004::fmd_32 fm(mod);
  return yosupo_193004::det_impl(a, fm);
 }
};

#endif // YOSUPO_193004_ENABLE
