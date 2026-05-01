#pragma once
#include "_common.hpp"
// =============================================================================
// Source: yosupo judge "Convolution (mod 998244353)" 最速提出を抽出。
//   提出: https://judge.yosupo.jp/submission/305538
// 抽出方針:
//   - I/O 部分 (Qinf/Qoutf, mmap など) は base.cpp に移譲して削除
//   - main() / auto_timer / std::mt19937_64 rng なども削除
//   - C++20 concepts (trivialT) は C++17 互換のため除去
//   - AVX2 intrinsics は simde 経由でも動くように切り替え
//   - 名前空間 yosupo_fastest_ntt で囲んで識別子衝突を避ける
//   - Conv::run(vector, vector) interface を末尾に追加
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#else
#pragma GCC target("avx2")
#include <immintrin.h>
#endif
#include <sys/mman.h>

// AVX-512 系ヘッダ由来のマクロ。clang + simde では未定義のことがあるので補完。
// _mm256_shuffle_epi32 の immediate (4 個の 2bit インデックスを 8bit に詰めたもの)。
#ifndef _MM_PERM_BBBB
#define _MM_PERM_BBBB 0x55  // (1,1,1,1)
#endif
#ifndef _MM_PERM_DDDD
#define _MM_PERM_DDDD 0xFF  // (3,3,3,3)
#endif

namespace yosupo_fastest_ntt {
using i32 = int32_t;
using i64 = int64_t;
using nu32 = uint32_t;
using nu64 = uint64_t;
using idt = std::size_t;
using I256 = __m256i;
template <class T, idt k> using arr = std::array<T, k>;

#define fn constexpr auto
#define let const auto
#define def inline auto
#define idef [[gnu::always_inline]] def

def store(void* p, I256 x) { _mm256_store_si256((I256*) p, x); }
def load256(const void* p) { return _mm256_load_si256((const I256*) p); }
def loadu256(const void* p) { return _mm256_loadu_si256((const I256*) p); }
fn shrk(nu32 x, nu32 M) { return std::min(x, x - M); }
fn dilt(nu32 x, nu32 M) { return std::min(x, x + M); }
fn reduce(nu64 x, nu32 ninv, nu32 M) -> nu32 { return (x + nu64(nu32(x) * ninv) * M) >> 32; }
fn mul(nu32 x, nu32 y, nu32 ninv, nu32 M) { return reduce(nu64(x) * y, ninv, M); }
fn mul_b_fixed(nu32 x, nu32 y, nu32 binv, nu32 M) -> nu32 { return (nu64(x) * y + nu64(binv * x) * M) >> 32; }
fn mul_s(nu32 x, nu32 y, nu32 ninv, nu32 M) { return shrk(reduce(nu64(x) * y, ninv, M), M); }
fn qpw(nu32 a, nu32 b, nu32 ninv, nu32 M, nu32 r) {
 for (; b; b >>= 1, a = mul(a, a, ninv, M)) {
  if (b & 1) { r = mul(r, a, ninv, M); }
 }
 return r;
}
fn qpw_s(nu32 a, nu32 b, nu32 ninv, nu32 M, nu32 r) { return shrk(qpw(a, b, ninv, M, r), M); }
def shrk32(I256 x, I256 M) { return _mm256_min_epu32(x, _mm256_sub_epi32(x, M)); }
def dilt32(I256 x, I256 M) { return _mm256_min_epu32(x, _mm256_add_epi32(x, M)); }
def add32(I256 x, I256 y) { return _mm256_add_epi32(x, y); }
def Lsub32(I256 x, I256 y, I256 M) { return _mm256_add_epi32(x, _mm256_sub_epi32(M, y)); }
def add32(I256 x, I256 y, I256 M) { return shrk32(_mm256_add_epi32(x, y), M); }
def sub32(I256 x, I256 y, I256 M) { return dilt32(_mm256_sub_epi32(x, y), M); }
def lmove(I256 x) { return _mm256_bsrli_epi128(x, 4); }
def reduce(I256 a, I256 b, I256 ninv, I256 M) {
 auto c = _mm256_mul_epu32(a, ninv), d = _mm256_mul_epu32(b, ninv);
 c = _mm256_mul_epu32(c, M), d = _mm256_mul_epu32(d, M);
 return _mm256_or_si256(lmove(_mm256_add_epi64(a, c)), _mm256_add_epi64(b, d));
}
template <int b_only_even = 0> def mul(I256 a, I256 b, I256 ninv, I256 M) {
 return reduce(_mm256_mul_epu32(a, b), _mm256_mul_epu32(lmove(a), b_only_even ? b : lmove(b)), ninv, M);
}
def mul_s(I256 a, I256 b, I256 ninv, I256 M) { return shrk32(mul(a, b, ninv, M), M); }
template <int b_only_even = 0> def mul_b_fixed(I256 a, I256 b, I256 bninv, I256 M) {
 I256 cc = _mm256_mul_epu32(a, bninv), dd = _mm256_mul_epu32(lmove(a), b_only_even ? bninv : lmove(bninv));
 I256 c = _mm256_mul_epu32(a, b), d = _mm256_mul_epu32(lmove(a), b_only_even ? b : lmove(b));
 cc = _mm256_mul_epu32(cc, M), dd = _mm256_mul_epu32(dd, M);
 return _mm256_or_si256(lmove(_mm256_add_epi64(c, cc)), _mm256_add_epi64(d, dd));
}
template <int b_only_even = 0> def mul_b_fixed_cross(I256 a, I256 b, I256 bninv, I256 M) {
 I256 cc = _mm256_mul_epu32(a, bninv), dd = _mm256_mul_epu32(lmove(a), b_only_even ? bninv : lmove(bninv));
 I256 c = _mm256_mul_epu32(a, b), d = _mm256_mul_epu32(lmove(a), b_only_even ? b : lmove(b));
 cc = _mm256_mul_epu32(cc, M), dd = _mm256_mul_epu32(dd, M);
 return _mm256_or_si256(_mm256_add_epi64(c, cc), lmove(_mm256_add_epi64(d, dd)));
}
def mul_upd_rt(I256 a, I256 bu, I256 M) {
 auto cc = _mm256_mul_epu32(a, bu), c = _mm256_mul_epu32(a, _mm256_srli_epi64(bu, 32));
 cc = _mm256_mul_epu32(cc, M);
 return shrk32(_mm256_srli_epi64(_mm256_add_epi64(c, cc), 32), M);
}
def mul_upd_rr(I256 a, I256 bu, I256 M) { return mul_b_fixed<1>(a, bu, _mm256_srli_epi64(bu, 32), M); }
idef expand(nu32 x) { return _mm256_set1_epi32(x); }
def select_load32(I256* f, idt p) { return ((nu32*) f)[p]; }
fn bcl(idt x) { return x < 2 ? 1 : idt(2) << (63 - __builtin_clzll(x - 1)); }

fn llimit = idt(1) << 20;
template <class T, idt al = 32> def _lalloc(idt n) { return new (std::align_val_t(al)) T[n]; }
template <class T, idt al = 32> def _lfree(T* p) { ::operator delete[](p, std::align_val_t(al)); }
template <class T, idt aln = 32> def lalloc(idt n) {
#if defined(__linux__) || defined(__unix__)
 idt bytes = sizeof(T) * n;
 if (bytes > llimit) {
  auto raw = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  madvise(raw, bytes, MADV_POPULATE_WRITE);
  return (T*) raw;
 }
#endif
 return _lalloc<T, aln>(n);
}
template <class T, idt aln = 32> def lfree(T* p, idt n) {
#if defined(__linux__) || defined(__unix__)
 idt bytes = sizeof(T) * n;
 if (bytes > llimit) {
  munmap(p, bytes);
  return;
 }
#endif
 _lfree<T, aln>(p);
}
template <class T> def ntt_copy(T* f, const T* g, idt n) { return (T*) std::memcpy(f, g, n * sizeof(T)); }
template <class T> def ntt_clear(T* f, idt n) { return (T*) std::memset(f, 0, n * sizeof(T)); }

namespace NTT_interal {
fn _lgmax = 26, _lg_iter_thresold = 6;
fn _iter_thresold = idt(1) << _lg_iter_thresold;
static_assert(_lg_iter_thresold % 2 == 0);
idef __conv8_4(nu32* __restrict__ f, nu32* __restrict__ g, std::array<nu32, 4> ww, I256 Ninv, I256 Mod, I256 Mod2) {
 alignas(64) nu32 awa[4][16];
 alignas(64) I256 res0[4] = {}, res1[4] = {};
#pragma GCC unroll(4)
 for (auto i = 0; i < 4; ++i) {
  store(g + i * 8, shrk32(shrk32(load256(g + i * 8), Mod2), Mod));
  auto ff = load256(f + i * 8);
  ff = shrk32(ff, Mod2);
  let ffw = shrk32(mul<1>(ff, expand(ww[i]), Ninv, Mod), Mod);
  ff = shrk32(ff, Mod);
  store(awa[i], ffw), store(awa[i] + 8, ff);
 }
 for (auto i = 0; i < 8; ++i) {
#pragma GCC unroll(4)
  for (auto j = 0; j < 4; ++j) {
   let bi = expand(g[j * 8 + i]);
   let aj = loadu256(awa[j] + 8 - i);
   let aj2 = lmove(aj);
   res0[j] = _mm256_add_epi64(res0[j], _mm256_mul_epu32(bi, aj));
   res1[j] = _mm256_add_epi64(res1[j], _mm256_mul_epu32(bi, aj2));
  }
 }
#pragma GCC unroll(4)
 for (auto i = 0; i < 4; ++i) { store(f + i * 8, shrk32(reduce(res0[i], res1[i], Ninv, Mod), Mod2)); }
}
struct NTT32_info {
 nu32 mod, mod2, ninv, one, r2, r3, imag, imagninv, RT1[_lgmax], RT3[_lgmax];
 alignas(32) arr<nu32, 8> rt3[_lgmax - 2], rt3i[_lgmax - 2], bwb, bwbi;
 nu64 rt4n[_lgmax - 3], rt4n2[_lgmax - 3], rt4ni[_lgmax - 3], rt4n2i[_lgmax - 3];
 alignas(32) arr<nu32, 8> rt4nr, rt4nr2, rt4nri, rt4nr2i;
 constexpr NTT32_info(nu32 m): mod(m), mod2(m * 2), ninv([&] { auto n= 2 + m; for(auto i= 0; i < 4; ++i) { n*= 2 + m * n; } return n; }()), one((-m) % m), r2((-nu64(m)) % m), r3(mul_s(r2, r2, ninv, m)), imag{}, imagninv{}, RT1{}, RT3{}, rt3{}, rt3i{}, bwb{}, bwbi{}, rt4n{}, rt4n2{}, rt4ni{}, rt4n2i{}, rt4nr{}, rt4nr2{}, rt4nri{}, rt4nr2i{} {
  auto k = __builtin_ctz(m - 1);
  auto _g = mul(3, r2, ninv, mod);
  for (;; ++_g) {
   if (qpw_s(_g, mod >> 1, ninv, mod, one) != one) { break; }
  }
  _g = qpw(_g, mod >> k, ninv, mod, one);
  nu32 rt1[_lgmax - 1] = {}, rt1i[_lgmax - 1] = {};
  rt1[k - 2] = _g, rt1i[k - 2] = qpw(_g, mod - 2, ninv, mod, one);
  for (auto i = k - 2; i > 0; --i) {
   rt1[i - 1] = mul(rt1[i], rt1[i], ninv, mod);
   rt1i[i - 1] = mul(rt1i[i], rt1i[i], ninv, mod);
  }
  RT1[k - 1] = qpw_s(_g, 3, ninv, mod, one);
  for (auto i = k - 1; i > 0; --i) { RT1[i - 1] = mul_s(RT1[i], RT1[i], ninv, mod); }
  imag = rt1[0], imagninv = imag * ninv;
  bwb = {rt1[1], 0, rt1[0], 0, mod - mul_s(rt1[0], rt1[1], ninv, mod)};
  bwbi = {rt1i[1], 0, rt1i[0], 0, mul_s(rt1i[0], rt1i[1], ninv, mod)};
  auto pr = one, pri = one;
  for (auto i = 0; i < k - 2; ++i) {
   let r = mul_s(pr, rt1[i + 1], ninv, mod), ri = mul_s(pri, rt1i[i + 1], ninv, mod);
   let r2 = mul_s(r, r, ninv, mod), r2i = mul_s(ri, ri, ninv, mod);
   let r3 = mul_s(r, r2, ninv, mod), r3i = mul_s(ri, r2i, ninv, mod);
   rt3[i] = {r * ninv, r, r2 * ninv, r2, r3 * ninv, r3};
   RT3[i + 2] = rt3[i][1];
   rt3i[i] = {ri * ninv, ri, r2i * ninv, r2i, r3i * ninv, r3i};
   pr = mul(pr, rt1i[i + 1], ninv, mod), pri = mul(pri, rt1[i + 1], ninv, mod);
  }
  nu32 w[8] = {}, wi[8] = {};
  w[0] = one, wi[0] = one;
  for (auto i = 0; i < 3; ++i) {
   pr = rt1[i], pri = rt1i[i];
   for (auto j = 1 << i, k = 0; k < j; ++k) {
    w[j + k] = mul_s(w[k], pr, ninv, mod);
    wi[j + k] = mul_s(wi[k], pri, ninv, mod);
   }
  }
  rt4nr = {w[7], w[0], w[6], w[1], w[5], w[2], w[4], w[3]};
  rt4nr2 = {w[1], w[3], w[1], w[0], w[0], w[2], w[0], w[1]};
  rt4nri = {wi[7], wi[0], wi[6], wi[1], wi[5], wi[2], wi[4], wi[3]};
  rt4nr2i = {wi[1], wi[3], wi[1], wi[0], wi[0], wi[2], wi[0], wi[1]};
  pr = one, pri = one;
  for (auto i = 1; i < k - 3; ++i) {
   let r = mul_s(pr, rt1[i + 2], ninv, mod), ri = mul_s(pri, rt1i[i + 2], ninv, mod);
   let r2 = mul_s(r, r, ninv, mod), r2i = mul_s(ri, ri, ninv, mod);
   let r4 = mul_s(r2, r2, ninv, mod), r4i = mul_s(r2i, r2i, ninv, mod);
   rt4n[i] = nu64(r * ninv) << 32 | r, rt4ni[i] = nu64(ri * ninv) << 32 | ri;
   rt4n2[i] = nu64(r2) << 32 | r4, rt4n2i[i] = nu64(r2i) << 32 | r4i;
   pr = mul(pr, rt1i[i + 2], ninv, mod), pri = mul(pri, rt1[i + 2], ninv, mod);
  }
 }
 def _vec_dif(I256* const f, idt n) const {
  alignas(32) arr<nu32, 8> st_1[_lgmax >> 1];
  let Mod = expand(mod), Mod2 = expand(mod2), Ninv = expand(ninv), Imag = expand(imag), ImagNinv = expand(imagninv);
  let id24 = _mm256_set_epi32(4, 0, 2, 0, 4, 0, 2, 0);
  let lgn = __builtin_ctzll(n);
  std::fill(st_1, st_1 + (lgn >> 1), bwb);
  let nn = n >> (lgn & 1), m = std::min(n, _iter_thresold);
  if (nn != n) {
   for (idt i = 0; i < nn; ++i) {
    auto p0 = f + i, p1 = f + nn + i;
    let f0 = load256(p0), f1 = load256(p1);
    let g0 = add32(f0, f1, Mod2), g1 = Lsub32(f0, f1, Mod2);
    store(p0, g0), store(p1, g1);
   }
  }
  for (auto L = nn >> 2; L > 0; L >>= 2) {
   for (idt i = 0; i < L; ++i) {
    auto p0 = f + i, p1 = p0 + L, p2 = p1 + L, p3 = p2 + L;
    let f1 = load256(p1), f3 = load256(p3), f2 = load256(p2), f0 = load256(p0);
    let g3 = mul_b_fixed<1>(Lsub32(f1, f3, Mod2), Imag, ImagNinv, Mod), g1 = add32(f1, f3, Mod2);
    let g0 = add32(f0, f2, Mod2), g2 = sub32(f0, f2, Mod2);
    let h0 = add32(g0, g1, Mod2), h1 = Lsub32(g0, g1, Mod2);
    let h2 = add32(g2, g3), h3 = Lsub32(g2, g3, Mod2);
    store(p0, h0), store(p1, h1), store(p2, h2), store(p3, h3);
   }
  }
  int t = std::min(_lg_iter_thresold, lgn) & -2, p = (t - 2) >> 1;
  for (idt j = 0; j < n; j += m, t = __builtin_ctzll(j) & -2, p = (t - 2) >> 1) {
   auto const g = f + j;
   for (idt l = (idt(1) << t), L = l >> 2; L > 1; l = L, L >>= 2, t -= 2, --p) {
    auto rt = load256(st_1 + p);
    for (idt i = (j == 0 ? l : 0), k = (j + i) >> t; i < m; i += l, ++k) {
     let r1 = _mm256_permutevar8x32_epi32(rt, id24);
     let r1Ninv = _mm256_permutevar8x32_epi32(_mm256_mul_epu32(rt, Ninv), id24);
     rt = mul_upd_rt(rt, load256(rt3 + __builtin_ctzll(~k)), Mod);
     let r2 = _mm256_shuffle_epi32(r1, _MM_PERM_BBBB), nr3 = _mm256_shuffle_epi32(r1, _MM_PERM_DDDD);
     let r2Ninv = _mm256_shuffle_epi32(r1Ninv, _MM_PERM_BBBB), nr3Ninv = _mm256_shuffle_epi32(r1Ninv, _MM_PERM_DDDD);
     for (idt j = 0; j < L; ++j) {
      auto p0 = g + i + j, p1 = p0 + L, p2 = p1 + L, p3 = p2 + L;
      let f1 = load256(p1), f3 = load256(p3), f2 = load256(p2), f0 = load256(p0);
      let g1 = mul_b_fixed<1>(f1, r1, r1Ninv, Mod), ng3 = mul_b_fixed<1>(f3, nr3, nr3Ninv, Mod);
      let g2 = mul_b_fixed<1>(f2, r2, r2Ninv, Mod), g0 = shrk32(f0, Mod2);
      let h3 = mul_b_fixed<1>(add32(g1, ng3), Imag, ImagNinv, Mod), h1 = sub32(g1, ng3, Mod2);
      let h0 = add32(g0, g2, Mod2), h2 = sub32(g0, g2, Mod2);
      let o0 = add32(h0, h1), o1 = Lsub32(h0, h1, Mod2);
      let o2 = add32(h2, h3), o3 = Lsub32(h2, h3, Mod2);
      store(p0, o0), store(p1, o1), store(p2, o2), store(p3, o3);
     }
    }
    store(st_1 + p, rt);
   }
   { // L == 1
    auto rt = load256(st_1);
    for (idt i = j + (j == 0) * 4; i < j + m; i += 4) {
     let r1 = _mm256_permutevar8x32_epi32(rt, id24);
     rt = mul_upd_rt(rt, load256(rt3 + __builtin_ctzll(~i >> 2)), Mod);
     auto p0 = f + i, p1 = p0 + 1, p2 = p0 + 2, p3 = p0 + 3;
     let f1 = load256(p1), f3 = load256(p3), f2 = load256(p2), f0 = load256(p0);
     let r2 = _mm256_shuffle_epi32(r1, _MM_PERM_BBBB), nr3 = _mm256_shuffle_epi32(r1, _MM_PERM_DDDD);
     let g1 = mul<1>(f1, r1, Ninv, Mod), ng3 = mul<1>(f3, nr3, Ninv, Mod);
     let g2 = mul<1>(f2, r2, Ninv, Mod), g0 = shrk32(f0, Mod2);
     let h3 = mul_b_fixed<1>(add32(g1, ng3), Imag, ImagNinv, Mod), h1 = sub32(g1, ng3, Mod2);
     let h0 = add32(g0, g2, Mod2), h2 = sub32(g0, g2, Mod2);
     let o0 = add32(h0, h1), o1 = Lsub32(h0, h1, Mod2);
     let o2 = add32(h2, h3), o3 = Lsub32(h2, h3, Mod2);
     store(p0, o0), store(p1, o1), store(p2, o2), store(p3, o3);
    }
    store(st_1, rt);
   }
  }
 }
 def _vec_dit(I256* const f, idt n) const {
  alignas(32) arr<nu32, 8> st_1[_lgmax >> 1];
  let Mod = expand(mod), Mod2 = expand(mod2), Ninv = expand(ninv), Imag = expand(imag), ImagNinv = expand(imagninv);
  let id24 = _mm256_set_epi32(4, 0, 2, 0, 4, 0, 2, 0);
  let lgn = __builtin_ctzll(n);
  std::fill(st_1, st_1 + (lgn >> 1), bwbi);
  let nn = n >> (lgn & 1), m = std::min(n, _iter_thresold);
  let fx = mul_s(mod - ((mod - 1) >> lgn), r3, ninv, mod);
  let Fx = expand(fx), FxNinv = expand(fx * ninv);
  store(st_1, Fx);
  for (idt j = 0; j < n; j += m) {
   auto tt = __builtin_ctzll(j + m), t = 4, p = 1;
   { // L == 1, append coefficient
    auto rt = load256(st_1);
    for (idt i = j; i < j + m; i += 4) {
     let r1 = _mm256_permutevar8x32_epi32(rt, id24);
     rt = mul_upd_rt(rt, load256(rt3i + __builtin_ctzll(~i >> 2)), Mod);
     auto const p0 = f + i, p1 = p0 + 1, p2 = p1 + 1, p3 = p2 + 1;
     let f2 = load256(p2), f3 = load256(p3), f0 = load256(p0), f1 = load256(p1);
     let r2 = _mm256_shuffle_epi32(r1, _MM_PERM_BBBB), r3 = _mm256_shuffle_epi32(r1, _MM_PERM_DDDD);
     let g3 = mul_b_fixed<1>(Lsub32(f3, f2, Mod2), Imag, ImagNinv, Mod), g2 = add32(f2, f3, Mod2);
     let g0 = add32(f0, f1, Mod2), g1 = sub32(f0, f1, Mod2);
     let h2 = Lsub32(g0, g2, Mod2), h3 = Lsub32(g1, g3, Mod2);
     let h0 = add32(g0, g2), h1 = add32(g1, g3);
     let o2 = mul<1>(h2, r2, Ninv, Mod), o0 = mul_b_fixed<1>(h0, Fx, FxNinv, Mod);
     let o1 = mul<1>(h1, r1, Ninv, Mod), o3 = mul<1>(h3, r3, Ninv, Mod);
     store(p0, o0), store(p1, o1), store(p2, o2), store(p3, o3);
    }
    store(st_1, rt);
   }
   for (idt l = 16, L = 4; t <= tt; L = l, l <<= 2, t += 2, ++p) {
    idt diff = j + m - std::max(l, m), i = 0;
    I256 rt = load256(st_1 + p), *const g = f + diff;
    if (diff == 0) {
     if (l == n) {
      for (idt i = 0; i < L; ++i) {
       auto const p0 = f + i, p1 = p0 + L, p2 = p1 + L, p3 = p2 + L;
       let f2 = load256(p2), f3 = load256(p3), f0 = load256(p0), f1 = load256(p1);
       let g3 = mul_b_fixed<1>(Lsub32(f3, f2, Mod2), Imag, ImagNinv, Mod), g2 = add32(f2, f3, Mod2);
       let g0 = add32(f0, f1, Mod2), g1 = sub32(f0, f1, Mod2);
       let h0 = add32(g0, g2, Mod2), h1 = add32(g1, g3, Mod2);
       let h2 = sub32(g0, g2, Mod2), h3 = sub32(g1, g3, Mod2);
       let o0 = shrk32(h0, Mod), o1 = shrk32(h1, Mod);
       let o2 = shrk32(h2, Mod), o3 = shrk32(h3, Mod);
       store(p0, o0), store(p1, o1), store(p2, o2), store(p3, o3);
      }
     } else {
      for (idt i = 0; i < L; ++i) {
       auto const p0 = f + i, p1 = p0 + L, p2 = p1 + L, p3 = p2 + L;
       let f2 = load256(p2), f3 = load256(p3), f0 = load256(p0), f1 = load256(p1);
       let g3 = mul_b_fixed<1>(Lsub32(f3, f2, Mod2), Imag, ImagNinv, Mod), g2 = add32(f2, f3, Mod2);
       let g0 = add32(f0, f1, Mod2), g1 = sub32(f0, f1, Mod2);
       let h0 = add32(g0, g2, Mod2), h1 = add32(g1, g3, Mod2);
       let h2 = sub32(g0, g2, Mod2), h3 = sub32(g1, g3, Mod2);
       store(p0, h0), store(p1, h1), store(p2, h2), store(p3, h3);
      }
     }
     i = l;
    }
    for (idt k = (j + i) >> t; i < m; i += l, ++k) {
     let r1 = _mm256_permutevar8x32_epi32(rt, id24);
     let r1Ninv = _mm256_permutevar8x32_epi32(_mm256_mul_epu32(rt, Ninv), id24);
     rt = mul_upd_rt(rt, load256(rt3i + __builtin_ctzll(~k)), Mod);
     let r2 = _mm256_shuffle_epi32(r1, _MM_PERM_BBBB), r3 = _mm256_shuffle_epi32(r1, _MM_PERM_DDDD);
     let r2Ninv = _mm256_shuffle_epi32(r1Ninv, _MM_PERM_BBBB), r3Ninv = _mm256_shuffle_epi32(r1Ninv, _MM_PERM_DDDD);
     for (idt j = 0; j < L; ++j) {
      auto const p0 = g + i + j, p1 = p0 + L, p2 = p1 + L, p3 = p2 + L;
      let f2 = load256(p2), f3 = load256(p3), f0 = load256(p0), f1 = load256(p1);
      let g3 = mul_b_fixed<1>(Lsub32(f3, f2, Mod2), Imag, ImagNinv, Mod), g2 = add32(f2, f3, Mod2);
      let g0 = add32(f0, f1, Mod2), g1 = sub32(f0, f1, Mod2);
      let h2 = Lsub32(g0, g2, Mod2), h3 = Lsub32(g1, g3, Mod2);
      let h0 = add32(g0, g2), h1 = add32(g1, g3);
      let o2 = mul_b_fixed<1>(h2, r2, r2Ninv, Mod), o0 = shrk32(h0, Mod2);
      let o1 = mul_b_fixed<1>(h1, r1, r1Ninv, Mod), o3 = mul_b_fixed<1>(h3, r3, r3Ninv, Mod);
      store(p0, o0), store(p1, o1), store(p2, o2), store(p3, o3);
     }
    }
    store(st_1 + p, rt);
   }
  }
  if (nn != n) {
   for (idt i = 0; i < nn; ++i) {
    auto const p0 = f + i, p1 = f + nn + i;
    let f0 = load256(p0), f1 = load256(p1);
    let g0 = add32(f0, f1, Mod2), g1 = sub32(f0, f1, Mod2);
    let h0 = shrk32(g0, Mod), h1 = shrk32(g1, Mod);
    store(p0, h0), store(p1, h1);
   }
  }
 }
 def _vec_cvdt8(I256* __restrict__ f, I256* __restrict__ g, idt lm) const {
  auto RR = one;
  let Ninv = expand(ninv), Mod = expand(mod), Mod2 = expand(mod2);
  for (idt i = 0; i < lm; i += 4) {
   let RRi = mul_b_fixed(RR, imag, imagninv, mod);
   __conv8_4((nu32*) (f + i), (nu32*) (g + i), {RR, mod2 - RR, RRi, mod2 - RRi}, Ninv, Mod, Mod2);
   RR = mul(RR, RT3[__builtin_ctzll(i + 4)], ninv, mod);
  }
 }
};
} // namespace NTT_interal

#undef fn
#undef let
#undef def
#undef idef

} // namespace yosupo_fastest_ntt

struct Conv {
 static vector<u32> run(const vector<u32>& a, const vector<u32>& b) {
  using namespace yosupo_fastest_ntt;
  using namespace yosupo_fastest_ntt::NTT_interal;
  static constexpr NTT32_info fntt(998244353);
  idt n = a.size(), m = b.size();
  if (!n || !m) return {};
  idt lm = bcl(std::max<idt>(64, n + m - 1));
  auto f = lalloc<nu32>(lm), g = lalloc<nu32>(lm);
  ntt_copy(f, a.data(), n);
  ntt_clear(f + n, lm - n);
  ntt_copy(g, b.data(), m);
  ntt_clear(g + m, lm - m);
  fntt._vec_dif((I256*) f, lm >> 3);
  fntt._vec_dif((I256*) g, lm >> 3);
  fntt._vec_cvdt8((I256*) f, (I256*) g, lm >> 3);
  fntt._vec_dit((I256*) f, lm >> 3);
  vector<u32> r(f, f + n + m - 1);
  lfree(f, lm);
  lfree(g, lm);
  return r;
 }
};
