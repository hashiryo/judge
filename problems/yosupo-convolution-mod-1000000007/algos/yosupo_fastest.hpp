#pragma once
#include "_common.hpp"
// =============================================================================
// Source: yosupo judge "Convolution (mod 10^9+7)" 最速提出を抽出。
//   提出: https://judge.yosupo.jp/submission/338813
//   元コードは cp-algorithms-aux (cp-algo) 由来:
//     https://github.com/cp-algorithms/cp-algorithms-aux
//   方式: 複素 FFT (double) + Karatsuba 3-split で arbitrary mod 畳み込み。
// 抽出方針:
//   - I/O 部分 (fastio mmap reader/writer) は base.cpp に移譲して削除
//   - main() 削除
//   - checkpoint デバッグ機能は no-op stub に
//   - cp_algo 名前空間はそのまま保持
//   - Conv::run interface を末尾に追加して std::vector<u32> ↔ cp_algo の base 変換
//   - C++23 (gnu++23) を前提とする
// =============================================================================

#pragma GCC optimize("Ofast,unroll-loops")
#include <bits/stdc++.h>
#include <ranges>
#include <bit>
// 元の cp_algo コードは `namespace stdx = std::experimental;` を作るが実際には未使用。
// gcc-14 + aarch64 で <experimental/simd> を include すると arm_neon.h の float16_t と
// <stdfloat> の std::float16_t が衝突するので、ここでは include しない。
// USE_SIMDE 環境では __AVX2__ が立たないので明示的に立てる (fallback 経路だと WA)。
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#ifndef __AVX2__
#define __AVX2__ 1
#endif
#elif defined(__AVX2__)
#include <immintrin.h>
#endif
#if defined(__linux__) || defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#define CP_ALGO_USE_MMAP 1
#include <sys/mman.h>
#else
#define CP_ALGO_USE_MMAP 0
#endif

namespace cp_algo {
// checkpoint: no-op (元実装はデバッグ用 timer)
template <bool final = false> inline void checkpoint(std::string const& = "") {}

namespace random {
inline uint64_t rng() {
 static std::mt19937_64 r(std::chrono::steady_clock::now().time_since_epoch().count());
 return r();
}
} // namespace random

template <typename T, std::size_t Align = 32> class big_alloc {
 static_assert(Align >= alignof(void*));
 static_assert(std::popcount(Align) == 1);
public:
 using value_type = T;
 template <class U> struct rebind { using other = big_alloc<U, Align>; };
 big_alloc() noexcept = default;
 template <typename U, std::size_t A> big_alloc(const big_alloc<U, A>&) noexcept {}
 [[nodiscard]] T* allocate(std::size_t n) {
  std::size_t padded = round_up(n * sizeof(T));
  std::size_t align = std::max<std::size_t>(alignof(T), Align);
#if CP_ALGO_USE_MMAP
  if (padded >= MEGABYTE) {
   void* raw = mmap(nullptr, padded, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#ifdef MADV_HUGEPAGE
   madvise(raw, padded, MADV_HUGEPAGE);
#endif
#ifdef MADV_POPULATE_WRITE
   madvise(raw, padded, MADV_POPULATE_WRITE);
#endif
   return static_cast<T*>(raw);
  }
#endif
  return static_cast<T*>(::operator new(padded, std::align_val_t(align)));
 }
 void deallocate(T* p, std::size_t n) noexcept {
  if (!p) return;
  std::size_t padded = round_up(n * sizeof(T));
  std::size_t align = std::max<std::size_t>(alignof(T), Align);
#if CP_ALGO_USE_MMAP
  if (padded >= MEGABYTE) { munmap(p, padded); return; }
#endif
  ::operator delete(p, padded, std::align_val_t(align));
 }
private:
 static constexpr std::size_t MEGABYTE = 1 << 20;
 static constexpr std::size_t round_up(std::size_t x) noexcept { return (x + Align - 1) / Align * Align; }
};

template <typename T, size_t len> using simd [[gnu::vector_size(len * sizeof(T))]] = T;
using i64x4 = simd<int64_t, 4>;
using u64x4 = simd<uint64_t, 4>;
using u32x8 = simd<uint32_t, 8>;
using i32x4 = simd<int32_t, 4>;
using u32x4 = simd<uint32_t, 4>;
using dx4 = simd<double, 4>;

[[gnu::always_inline]] inline dx4 abs(dx4 a) { return a < 0 ? -a : a; }

static constexpr dx4 magic_pad = dx4() + (3ULL << 51);
[[gnu::always_inline]] inline i64x4 lround(dx4 x) { return i64x4(x + magic_pad) - i64x4(magic_pad); }
[[gnu::always_inline]] inline dx4 to_double(i64x4 x) { return dx4(x + i64x4(magic_pad)) - magic_pad; }

[[gnu::always_inline]] inline dx4 round(dx4 a) {
 return dx4{std::nearbyint(a[0]), std::nearbyint(a[1]), std::nearbyint(a[2]), std::nearbyint(a[3])};
}

[[gnu::always_inline]] inline u64x4 montgomery_reduce(u64x4 x, u64x4 mod, u64x4 imod) {
 auto x_ninv = u64x4(u32x8(x) * u32x8(imod));
#ifdef __AVX2__
 x += u64x4(_mm256_mul_epu32(__m256i(x_ninv), __m256i(mod)));
#else
 x += x_ninv * mod;
#endif
 return x >> 32;
}
[[gnu::always_inline]] inline u64x4 montgomery_mul(u64x4 x, u64x4 y, u64x4 mod, u64x4 imod) {
#ifdef __AVX2__
 return montgomery_reduce(u64x4(_mm256_mul_epu32(__m256i(x), __m256i(y))), mod, imod);
#else
 return montgomery_reduce(x * y, mod, imod);
#endif
}
[[gnu::always_inline]] inline dx4 rotate_right(dx4 x) {
#if defined(__clang__)
 return __builtin_shufflevector(x, x, 3, 0, 1, 2);
#else
 static constexpr u64x4 shuffler = {3, 0, 1, 2};
 return __builtin_shuffle(x, shuffler);
#endif
}

template <typename T> struct complex {
 using value_type = T;
 T x, y;
 constexpr complex(): x(), y() {}
 constexpr complex(T x): x(x), y() {}
 constexpr complex(T x, T y): x(x), y(y) {}
 complex& operator*=(T t) { x*= t; y*= t; return *this; }
 complex& operator/=(T t) { x/= t; y/= t; return *this; }
 complex operator*(T t) const { return complex(*this)*= t; }
 complex operator/(T t) const { return complex(*this)/= t; }
 complex& operator+=(complex t) { x+= t.x; y+= t.y; return *this; }
 complex& operator-=(complex t) { x-= t.x; y-= t.y; return *this; }
 complex operator*(complex t) const { return {x * t.x - y * t.y, x * t.y + y * t.x}; }
 complex operator/(complex t) const { return *this * t.conj() / t.norm(); }
 complex operator+(complex t) const { return complex(*this)+= t; }
 complex operator-(complex t) const { return complex(*this)-= t; }
 complex& operator*=(complex t) { return *this= *this * t; }
 complex& operator/=(complex t) { return *this= *this / t; }
 complex operator-() const { return {-x, -y}; }
 complex conj() const { return {x, -y}; }
 T norm() const { return x * x + y * y; }
 T abs() const { return std::sqrt(norm()); }
 T const real() const { return x; }
 T const imag() const { return y; }
 T& real() { return x; }
 T& imag() { return y; }
 static constexpr complex polar(T r, T theta) { return {T(r * cos(theta)), T(r * sin(theta))}; }
 auto operator<=>(complex const&) const = default;
};
template <typename T> complex<T> operator*(auto x, complex<T> y) { return y*= x; }
template <typename T> complex<T> conj(complex<T> x) { return x.conj(); }
template <typename T> T norm(complex<T> x) { return x.norm(); }
template <typename T> T abs(complex<T> x) { return x.abs(); }
template <typename T> T& real(complex<T>& x) { return x.real(); }
template <typename T> T& imag(complex<T>& x) { return x.imag(); }
template <typename T> T const real(complex<T> const& x) { return x.real(); }
template <typename T> T const imag(complex<T> const& x) { return x.imag(); }
template <typename T> constexpr complex<T> polar(T r, T theta) { return complex<T>::polar(r, theta); }

namespace math {
const int maxn = 1 << 19;
const int magic = 64;
auto bpow(auto const& x, auto n, auto const& one, auto op) {
 if (n == 0) return one;
 auto t = bpow(x, n / 2, one, op);
 t = op(t, t);
 if (n % 2) t = op(t, x);
 return t;
}
auto bpow(auto x, auto n, auto ans) { return bpow(x, n, ans, std::multiplies{}); }
template <typename T> T bpow(T const& x, auto n) { return bpow(x, n, T(1)); }

template <typename modint, typename _Int> struct modint_base {
 using Int = _Int;
 using UInt = std::make_unsigned_t<Int>;
 static constexpr size_t bits = sizeof(Int) * 8;
 using Int2 = std::conditional_t<bits <= 32, int64_t, __int128_t>;
 using UInt2 = std::conditional_t<bits <= 32, uint64_t, __uint128_t>;
 static Int mod() { return modint::mod(); }
 static Int remod() { return modint::remod(); }
 static UInt2 modmod() { return UInt2(mod()) * mod(); }
 modint_base(): r(0) {}
 modint_base(Int2 rr) { to_modint().setr(UInt((rr + modmod()) % mod())); }
 modint inv() const { return bpow(to_modint(), mod() - 2); }
 modint operator-() const { modint neg; neg.r = std::min(-r, remod() - r); return neg; }
 modint& operator/=(const modint& t) { return to_modint()*= t.inv(); }
 modint& operator*=(const modint& t) { r = UInt(UInt2(r) * t.r % mod()); return to_modint(); }
 modint& operator+=(const modint& t) { r+= t.r; r = std::min(r, r - remod()); return to_modint(); }
 modint& operator-=(const modint& t) { r-= t.r; r = std::min(r, r + remod()); return to_modint(); }
 modint operator+(const modint& t) const { return modint(to_modint())+= t; }
 modint operator-(const modint& t) const { return modint(to_modint())-= t; }
 modint operator*(const modint& t) const { return modint(to_modint())*= t; }
 modint operator/(const modint& t) const { return modint(to_modint())/= t; }
 auto operator==(const modint& t) const { return to_modint().getr() == t.getr(); }
 auto operator!=(const modint& t) const { return to_modint().getr() != t.getr(); }
 Int rem() const { UInt R = to_modint().getr(); return R - (R > (UInt) mod() / 2) * mod(); }
 void setr(UInt rr) { r = rr; }
 UInt getr() const { return r; }
 static UInt modmod8() { return UInt(8 * modmod()); }
 void add_unsafe(UInt t) { r+= t; }
 void pseudonormalize() { r = std::min(r, r - modmod8()); }
 modint const& normalize() {
  if (r >= (UInt) mod()) r%= mod();
  return to_modint();
 }
 void setr_direct(UInt rr) { r = rr; }
 UInt getr_direct() const { return r; }
protected:
 UInt r;
private:
 modint& to_modint() { return static_cast<modint&>(*this); }
 modint const& to_modint() const { return static_cast<modint const&>(*this); }
};
template <typename modint> concept modint_type = std::is_base_of_v<modint_base<modint, typename modint::Int>, modint>;
template <auto m> struct modint: modint_base<modint<m>, decltype(m)> {
 using Base = modint_base<modint<m>, decltype(m)>;
 using Base::Base;
 static constexpr typename Base::Int mod() { return m; }
 static constexpr typename Base::UInt remod() { return m; }
 auto getr() const { return Base::r; }
};
inline constexpr auto inv2(auto x) {
 std::make_unsigned_t<decltype(x)> y = 1;
 while (y * x != 1) y*= 2 - x * y;
 return y;
}
} // namespace math
} // namespace cp_algo

namespace cp_algo::math::fft {
namespace stdx_dummy {} // silence unused
static constexpr size_t flen = 4;
using ftype = double;
using vftype = dx4;
using point = complex<ftype>;
using vpoint = complex<vftype>;
static constexpr vftype vz = {};
inline vpoint vi(vpoint const& r) { return {-imag(r), real(r)}; }

struct cvector {
 std::vector<vpoint, big_alloc<vpoint>> r;
 cvector(size_t n) {
  n = std::max(flen, std::bit_ceil(n));
  r.resize(n / flen);
  checkpoint("cvector create");
 }
 vpoint& at(size_t k) { return r[k / flen]; }
 vpoint at(size_t k) const { return r[k / flen]; }
 template <class pt = point> void set(size_t k, pt t) {
  if constexpr (std::is_same_v<pt, point>) {
   real(r[k / flen])[k % flen] = real(t);
   imag(r[k / flen])[k % flen] = imag(t);
  } else { at(k) = t; }
 }
 template <class pt = point> pt get(size_t k) const {
  if constexpr (std::is_same_v<pt, point>) {
   return {real(r[k / flen])[k % flen], imag(r[k / flen])[k % flen]};
  } else { return at(k); }
 }
 size_t size() const { return flen * r.size(); }
 static constexpr size_t eval_arg(size_t n) {
  if (n < pre_evals) return eval_args[n];
  return eval_arg(n / 2) | (n & 1) << (std::bit_width(n) - 1);
 }
 // libc++ では cos/sin が非 constexpr なので constexpr を外す (libstdc++ では問題ない)。
 static point eval_point(size_t n) {
  if (n % 2) return -eval_point(n - 1);
  if (n % 4) return eval_point(n - 2) * point(0, 1);
  if (n / 4 < pre_evals) return evalp[n / 4];
  return polar<ftype>(1., std::numbers::pi / (ftype) std::bit_floor(n) * (ftype) eval_arg(n));
 }
 inline static const std::array<point, 32> roots = []() {
  std::array<point, 32> res{};
  for (size_t i = 2; i < 32; i++) res[i] = polar<ftype>(1., std::numbers::pi / (1ull << (i - 2)));
  return res;
 }();
 static point root(size_t n) { return roots[std::bit_width(n)]; }
 template <int step> static void exec_on_eval(size_t n, size_t k, auto&& callback) {
  callback(k, root(4 * step * n) * eval_point(step * k));
 }
 template <int step> static void exec_on_evals(size_t n, auto&& callback) {
  point factor = root(4 * step * n);
  for (size_t i = 0; i < n; i++) callback(i, factor * eval_point(step * i));
 }
 void dot(cvector const& t) {
  size_t n = this->size();
  exec_on_evals<1>(n / flen, [&](size_t k, point rt) {
   k*= flen;
   auto [Ax, Ay] = at(k);
   auto Bv = t.at(k);
   vpoint res = vz;
   for (size_t i = 0; i < flen; i++) {
    res+= vpoint(vz + Ax[i], vz + Ay[i]) * Bv;
    real(Bv) = rotate_right(real(Bv));
    imag(Bv) = rotate_right(imag(Bv));
    auto x_ = real(Bv)[0], y_ = imag(Bv)[0];
    real(Bv)[0] = x_ * real(rt) - y_ * imag(rt);
    imag(Bv)[0] = x_ * imag(rt) + y_ * real(rt);
   }
   set(k, res);
  });
  checkpoint("dot");
 }
 void ifft() {
  const size_t n = size();
  [[gnu::assume(n <= 1 << 20)]];
  bool parity = std::countr_zero(n) % 2;
  if (parity) {
   exec_on_evals<2>(n / (2 * flen), [&](size_t k, point rt) {
    k*= 2 * flen;
    vpoint cvrt = {vz + real(rt), vz - imag(rt)};
    auto B = at(k) - at(k + flen);
    at(k)+= at(k + flen);
    at(k + flen) = B * cvrt;
   });
  }
  for (size_t leaf = 3 * flen; leaf < n; leaf+= 4 * flen) {
   size_t level = std::countr_one(leaf + 3);
   for (size_t lvl = 4 + parity; lvl <= level; lvl+= 2) {
    size_t i = (1 << lvl) / 4;
    exec_on_eval<4>(n >> lvl, leaf >> lvl, [&](size_t k, point rt) {
     k <<= lvl;
     vpoint v1 = {vz + real(rt), vz - imag(rt)};
     vpoint v2 = v1 * v1;
     vpoint v3 = v1 * v2;
     for (size_t j = k; j < k + i; j+= flen) {
      auto A = at(j);
      auto B = at(j + i);
      auto C = at(j + 2 * i);
      auto D = at(j + 3 * i);
      at(j) = ((A + B) + (C + D));
      at(j + 2 * i) = ((A + B) - (C + D)) * v2;
      at(j + i) = ((A - B) - vi(C - D)) * v1;
      at(j + 3 * i) = ((A - B) + vi(C - D)) * v3;
     }
    });
   }
  }
  checkpoint("ifft");
  for (size_t k = 0; k < n; k+= flen) set(k, get<vpoint>(k)/= vz + (ftype)(n / flen));
 }
 void fft() {
  const size_t n = size();
  [[gnu::assume(n <= 1 << 20)]];
  bool parity = std::countr_zero(n) % 2;
  for (size_t leaf = 0; leaf < n; leaf+= 4 * flen) {
   size_t level = std::countr_zero(n + leaf);
   level-= level % 2 != parity;
   for (size_t lvl = level; lvl >= 4; lvl-= 2) {
    size_t i = (1 << lvl) / 4;
    exec_on_eval<4>(n >> lvl, leaf >> lvl, [&](size_t k, point rt) {
     k <<= lvl;
     vpoint v1 = {vz + real(rt), vz + imag(rt)};
     vpoint v2 = v1 * v1;
     vpoint v3 = v1 * v2;
     for (size_t j = k; j < k + i; j+= flen) {
      auto A = at(j);
      auto B = at(j + i) * v1;
      auto C = at(j + 2 * i) * v2;
      auto D = at(j + 3 * i) * v3;
      at(j) = (A + C) + (B + D);
      at(j + i) = (A + C) - (B + D);
      at(j + 2 * i) = (A - C) + vi(B - D);
      at(j + 3 * i) = (A - C) - vi(B - D);
     }
    });
   }
  }
  if (parity) {
   exec_on_evals<2>(n / (2 * flen), [&](size_t k, point rt) {
    k*= 2 * flen;
    vpoint vrt = {vz + real(rt), vz + imag(rt)};
    auto t = at(k + flen) * vrt;
    at(k + flen) = at(k) - t;
    at(k)+= t;
   });
  }
  checkpoint("fft");
 }
 static constexpr size_t pre_evals = 1 << 16;
 static const std::array<size_t, pre_evals> eval_args;
 static const std::array<point, pre_evals> evalp;
};
inline const std::array<size_t, cvector::pre_evals> cvector::eval_args = []() {
 std::array<size_t, cvector::pre_evals> res = {};
 for (size_t i = 1; i < cvector::pre_evals; i++) res[i] = res[i >> 1] | (i & 1) << (std::bit_width(i) - 1);
 return res;
}();
inline const std::array<point, cvector::pre_evals> cvector::evalp = []() {
 std::array<point, cvector::pre_evals> res = {};
 res[0] = 1;
 for (size_t n = 1; n < cvector::pre_evals; n++)
  res[n] = polar<ftype>(1., std::numbers::pi * ftype(cvector::eval_args[n]) / ftype(4 * std::bit_floor(n)));
 return res;
}();

template <modint_type base> struct dft {
 cvector A, B;
 static base factor, ifactor;
 using Int2 = typename base::Int2;
 static bool _init;
 static int split() {
  static const int splt = int(std::sqrt(base::mod())) + 1;
  return splt;
 }
 static u64x4 mod_v, imod_v;
 static void init() {
  if (!_init) {
   factor = 1 + random::rng() % (base::mod() - 1);
   ifactor = base(1) / factor;
   mod_v = u64x4() + base::mod();
   imod_v = u64x4() + inv2(-base::mod());
   _init = true;
  }
 }
 dft(auto const& a, size_t n): A(n), B(n) {
  init();
  base b2x32 = bpow(base(2), 32);
  u64x4 cur = {(bpow(factor, 1) * b2x32).getr(), (bpow(factor, 2) * b2x32).getr(), (bpow(factor, 3) * b2x32).getr(), (bpow(factor, 4) * b2x32).getr()};
  u64x4 step4 = u64x4{} + (bpow(factor, 4) * b2x32).getr();
  u64x4 stepn = u64x4{} + (bpow(factor, n) * b2x32).getr();
  for (size_t i = 0; i < std::min(n, std::size(a)); i+= flen) {
   auto splt = [&](size_t i, auto mul) {
    if (i >= std::size(a)) return std::pair{vftype(), vftype()};
    u64x4 au = {i < std::size(a) ? a[i].getr() : 0, i + 1 < std::size(a) ? a[i + 1].getr() : 0, i + 2 < std::size(a) ? a[i + 2].getr() : 0, i + 3 < std::size(a) ? a[i + 3].getr() : 0};
    au = montgomery_mul(au, mul, mod_v, imod_v);
    au = au >= base::mod() ? au - base::mod() : au;
    auto ai = to_double(i64x4(au >= base::mod() / 2 ? au - base::mod() : au));
    auto quo = round(ai / split());
    return std::pair{ai - quo * split(), quo};
   };
   auto [rai, qai] = splt(i, cur);
   auto [rani, qani] = splt(n + i, montgomery_mul(cur, stepn, mod_v, imod_v));
   A.at(i) = vpoint(rai, rani);
   B.at(i) = vpoint(qai, qani);
   cur = montgomery_mul(cur, step4, mod_v, imod_v);
  }
  checkpoint("dft init");
  if (n) {
   A.fft();
   B.fft();
  }
 }
 void dot(auto&& C, auto const& D) {
  cvector::exec_on_evals<1>(A.size() / flen, [&](size_t k, point rt) {
   k*= flen;
   auto [Ax, Ay] = A.at(k);
   auto [Bx, By] = B.at(k);
   vpoint AC, AD, BC, BD;
   AC = AD = BC = BD = vz;
   auto Cv = C.at(k), Dv = D.at(k);
   for (size_t i = 0; i < flen; i++) {
    vpoint Av = {vz + Ax[i], vz + Ay[i]}, Bv = {vz + Bx[i], vz + By[i]};
    AC+= Av * Cv;
    AD+= Av * Dv;
    BC+= Bv * Cv;
    BD+= Bv * Dv;
    real(Cv) = rotate_right(real(Cv));
    imag(Cv) = rotate_right(imag(Cv));
    real(Dv) = rotate_right(real(Dv));
    imag(Dv) = rotate_right(imag(Dv));
    auto cx = real(Cv)[0], cy = imag(Cv)[0];
    auto dx_ = real(Dv)[0], dy = imag(Dv)[0];
    real(Cv)[0] = cx * real(rt) - cy * imag(rt);
    imag(Cv)[0] = cx * imag(rt) + cy * real(rt);
    real(Dv)[0] = dx_ * real(rt) - dy * imag(rt);
    imag(Dv)[0] = dx_ * imag(rt) + dy * real(rt);
   }
   A.at(k) = AC;
   C.at(k) = AD + BC;
   B.at(k) = BD;
  });
  checkpoint("dot");
 }
 void recover_mod(auto&& C, auto& res, size_t k) {
  size_t check = (k + flen - 1) / flen * flen;
  assert(res.size() >= check);
  size_t n = A.size();
  auto const splitsplit = base(split() * split()).getr();
  base b2x32 = bpow(base(2), 32);
  base b2x64 = bpow(base(2), 64);
  u64x4 cur = {(bpow(ifactor, 2) * b2x64).getr(), (bpow(ifactor, 3) * b2x64).getr(), (bpow(ifactor, 4) * b2x64).getr(), (bpow(ifactor, 5) * b2x64).getr()};
  u64x4 step4 = u64x4{} + (bpow(ifactor, 4) * b2x32).getr();
  u64x4 stepn = u64x4{} + (bpow(ifactor, n) * b2x32).getr();
  for (size_t i = 0; i < std::min(n, k); i+= flen) {
   auto [Ax, Ay] = A.at(i);
   auto [Bx, By] = B.at(i);
   auto [Cx, Cy] = C.at(i);
   auto set_i = [&](size_t i, auto A_, auto B_, auto C_, auto mul) {
    auto A0 = lround(A_), A1 = lround(C_), A2 = lround(B_);
    auto Ai = A0 + A1 * (int64_t)split() + A2 * (int64_t)splitsplit + (int64_t)(uint64_t)base::modmod();
    auto Au = montgomery_reduce(u64x4(Ai), mod_v, imod_v);
    Au = montgomery_mul(Au, mul, mod_v, imod_v);
    Au = Au >= base::mod() ? Au - base::mod() : Au;
    for (size_t j = 0; j < flen; j++) res[i + j].setr(typename base::UInt(Au[j]));
   };
   set_i(i, Ax, Bx, Cx, cur);
   if (i + n < k) set_i(i + n, Ay, By, Cy, montgomery_mul(cur, stepn, mod_v, imod_v));
   cur = montgomery_mul(cur, step4, mod_v, imod_v);
  }
  checkpoint("recover mod");
 }
 void mul(auto&& C, auto const& D, auto& res, size_t k) {
  assert(A.size() == C.size());
  size_t n = A.size();
  if (!n) {
   res = {};
   return;
  }
  dot(C, D);
  A.ifft();
  B.ifft();
  C.ifft();
  recover_mod(C, res, k);
 }
 void mul_inplace(auto&& B_, auto& res, size_t k) { mul(B_.A, B_.B, res, k); }
 void mul(auto const& B_, auto& res, size_t k) { mul(cvector(B_.A), B_.B, res, k); }
};
template <modint_type base> base dft<base>::factor = 1;
template <modint_type base> base dft<base>::ifactor = 1;
template <modint_type base> bool dft<base>::_init = false;
template <modint_type base> u64x4 dft<base>::mod_v = {};
template <modint_type base> u64x4 dft<base>::imod_v = {};

inline void mul_slow(auto& a, auto const& b, size_t k) {
 if (std::empty(a) || std::empty(b)) {
  a.clear();
 } else {
  size_t n = std::min(k, std::size(a));
  size_t m = std::min(k, std::size(b));
  a.resize(k);
  for (int j = int(k - 1); j >= 0; j--) {
   a[j]*= b[0];
   for (int i = std::max(j - (int) n, 0) + 1; i < std::min(j + 1, (int) m); i++) a[j]+= a[j - i] * b[i];
  }
 }
}
inline size_t com_size(size_t as, size_t bs) {
 if (!as || !bs) return 0;
 return std::max(flen, std::bit_ceil(as + bs - 1) / 2);
}
inline void mul_truncate(auto& a, auto const& b, size_t k) {
 using base = std::decay_t<decltype(a[0])>;
 if (std::min({k, std::size(a), std::size(b)}) < magic) {
  mul_slow(a, b, k);
  return;
 }
 auto n = std::max(flen, std::bit_ceil(std::min(k, std::size(a)) + std::min(k, std::size(b)) - 1) / 2);
 auto A = dft<base>(a | std::views::take(k), n);
 auto B = dft<base>(b | std::views::take(k), n);
 a.resize((k + flen - 1) / flen * flen);
 A.mul_inplace(B, a, k);
 a.resize(k);
}
inline void mod_split(auto&& x, size_t n, auto k) {
 using base = std::decay_t<decltype(k)>;
 dft<base>::init();
 assert(std::size(x) == 2 * n);
 u64x4 cur = u64x4{} + (k * bpow(base(2), 32)).getr();
 for (size_t i = 0; i < n; i+= flen) {
  u64x4 xl = *reinterpret_cast<u64x4*>(&x[i]);
  u64x4 xr = *reinterpret_cast<u64x4*>(&x[n + i]);
  xr = montgomery_mul(xr, cur, dft<base>::mod_v, dft<base>::imod_v);
  xr = xr >= base::mod() ? xr - base::mod() : xr;
  auto t = xr;
  xr = xl - t;
  xl+= t;
  xl = xl >= base::mod() ? xl - base::mod() : xl;
  xr = xr >= base::mod() ? xr + base::mod() : xr;
  for (size_t kk = 0; kk < flen; kk++) {
   x[i + kk].setr(typename base::UInt(xl[kk]));
   x[n + i + kk].setr(typename base::UInt(xr[kk]));
  }
 }
 cp_algo::checkpoint("mod split");
}
inline void cyclic_mul(auto& a, auto&& b, size_t k) {
 assert(std::popcount(k) == 1);
 assert(std::size(a) == std::size(b) && std::size(a) == k);
 using base = std::decay_t<decltype(a[0])>;
 dft<base>::init();
 if (k <= (1 << 16)) {
  auto ap = std::ranges::to<std::vector<base, big_alloc<base>>>(a);
  mul_truncate(ap, b, 2 * k);
  mod_split(ap, k, bpow(dft<base>::factor, k));
  std::ranges::copy(ap | std::views::take(k), begin(a));
  return;
 }
 k/= 2;
 auto factor = bpow(dft<base>::factor, k);
 mod_split(a, k, factor);
 mod_split(b, k, factor);
 auto la = std::span(a).first(k);
 auto lb = std::span(b).first(k);
 auto ra = std::span(a).last(k);
 auto rb = std::span(b).last(k);
 cyclic_mul(la, lb, k);
 auto A = dft<base>(ra, k / 2);
 auto B = dft<base>(rb, k / 2);
 A.mul_inplace(B, ra, k);
 base i2 = base(2).inv();
 factor = factor.inv() * i2;
 for (size_t i = 0; i < k; i++) {
  auto t = (a[i] + a[i + k]) * i2;
  a[i + k] = (a[i] - a[i + k]) * factor;
  a[i] = t;
 }
 cp_algo::checkpoint("mod join");
}
inline auto make_copy(auto&& x) { return x; }
inline void cyclic_mul(auto& a, auto const& b, size_t k) { return cyclic_mul(a, make_copy(b), k); }
inline void mul(auto& a, auto&& b) {
 size_t N = size(a) + size(b);
 if (N > (1 << 20)) {
  N--;
  size_t NN = std::bit_ceil(N);
  a.resize(NN);
  b.resize(NN);
  cyclic_mul(a, b, NN);
  a.resize(N);
 } else {
  mul_truncate(a, b, N - 1);
 }
}
inline void mul(auto& a, auto const& b) {
 size_t N = size(a) + size(b);
 if (N > (1 << 20)) {
  mul(a, make_copy(b));
 } else {
  mul_truncate(a, b, N - 1);
 }
}
} // namespace cp_algo::math::fft

struct Conv {
 using base = cp_algo::math::modint<(int64_t)1'000'000'007>;
 static vector<u32> run(const vector<u32>& a, const vector<u32>& b) {
  if (a.empty() || b.empty()) return {};
  std::vector<base, cp_algo::big_alloc<base>> aa(a.size()), bb(b.size());
  for (size_t i = 0; i < a.size(); ++i) aa[i].setr(a[i]);
  for (size_t i = 0; i < b.size(); ++i) bb[i].setr(b[i]);
  cp_algo::math::fft::mul(aa, bb);
  vector<u32> r(aa.size());
  for (size_t i = 0; i < aa.size(); ++i) r[i] = (u32) aa[i].getr();
  return r;
 }
};
