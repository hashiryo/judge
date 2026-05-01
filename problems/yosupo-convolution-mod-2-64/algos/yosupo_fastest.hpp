#pragma once
#include "_common.hpp"
// =============================================================================
// Source: yosupo judge "Convolution (mod 2^64)" 最速提出を抽出。
//   提出: https://judge.yosupo.jp/submission/337130
//   元コードは cp-algorithms-aux (cp-algo) 由来。
//   方式: 複素 FFT (double) + u64 を 4 × int16 に分割して畳み込み。
//   parts=4 系列を全部 FFT/IFFT し、各 chunk を u64 ラップで合成して結果。
// 抽出方針:
//   - I/O 部分 (blazingio, fastio) は base.cpp に移譲して削除
//   - main() / solve() 削除
//   - checkpoint デバッグは no-op stub に
//   - cp_algo 名前空間はそのまま保持
//   - <generator> 関連 (big_generator, ranges deduction guide) は使わないので削除
//   - Conv::run interface を末尾に追加
//   - C++23 (gnu++23) を前提とする
// =============================================================================

#pragma GCC optimize("Ofast,unroll-loops")
#include <bits/stdc++.h>
#include <ranges>
#include <bit>
#if __has_include(<experimental/simd>)
#include <experimental/simd>
#endif

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
template <bool final = false> inline void checkpoint(auto const&) {}
template <bool final = false> inline void checkpoint() {}

namespace random {
inline uint64_t rng() {
 static std::mt19937_64 g(std::chrono::steady_clock::now().time_since_epoch().count());
 return g();
}
} // namespace random

template <typename T, std::size_t Align = 32> class big_alloc {
 static_assert(Align >= alignof(void*));
 static_assert(std::popcount(Align) == 1);
public:
 using value_type = T;
 template <class U> struct rebind { using other = big_alloc<U, Align>; };
 constexpr bool operator==(const big_alloc&) const = default;
 constexpr bool operator!=(const big_alloc&) const = default;
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
template <typename T> using big_vector = std::vector<T, big_alloc<T>>;

template <typename T, size_t len> using simd [[gnu::vector_size(len * sizeof(T))]] = T;
using i64x4 = simd<int64_t, 4>;
using u64x4 = simd<uint64_t, 4>;
using u32x8 = simd<uint32_t, 8>;
using i32x4 = simd<int32_t, 4>;
using u32x4 = simd<uint32_t, 4>;
using i16x4 = simd<int16_t, 4>;
using u8x32 = simd<uint8_t, 32>;
using dx4 = simd<double, 4>;

inline dx4 abs(dx4 a) { return dx4{std::abs(a[0]), std::abs(a[1]), std::abs(a[2]), std::abs(a[3])}; }
static constexpr dx4 magic_pad = dx4() + (3ULL << 51);
inline i64x4 lround(dx4 x) { return i64x4(x + magic_pad) - i64x4(magic_pad); }
inline dx4 to_double(i64x4 x) { return dx4(x + i64x4(magic_pad)) - magic_pad; }
inline dx4 round(dx4 a) {
 return dx4{std::nearbyint(a[0]), std::nearbyint(a[1]), std::nearbyint(a[2]), std::nearbyint(a[3])};
}
inline dx4 rotate_right(dx4 x) {
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
 inline constexpr complex(): x(), y() {}
 inline constexpr complex(T const& x): x(x), y() {}
 inline constexpr complex(T const& x, T const& y): x(x), y(y) {}
 inline complex& operator*=(T const& t) { x*= t; y*= t; return *this; }
 inline complex& operator/=(T const& t) { x/= t; y/= t; return *this; }
 inline complex operator*(T const& t) const { return complex(*this)*= t; }
 inline complex operator/(T const& t) const { return complex(*this)/= t; }
 inline complex& operator+=(complex const& t) { x+= t.x; y+= t.y; return *this; }
 inline complex& operator-=(complex const& t) { x-= t.x; y-= t.y; return *this; }
 inline complex operator*(complex const& t) const { return {x * t.x - y * t.y, x * t.y + y * t.x}; }
 inline complex operator+(complex const& t) const { return complex(*this)+= t; }
 inline complex operator-(complex const& t) const { return complex(*this)-= t; }
 inline complex& operator*=(complex const& t) { return *this= *this * t; }
 inline complex operator-() const { return {-x, -y}; }
 inline complex conj() const { return {x, -y}; }
 inline T norm() const { return x * x + y * y; }
 inline T abs() const { return std::sqrt(norm()); }
 inline T const real() const { return x; }
 inline T const imag() const { return y; }
 inline T& real() { return x; }
 inline T& imag() { return y; }
 inline static constexpr complex polar(T r, T theta) { return {T(r * cos(theta)), T(r * sin(theta))}; }
 inline auto operator<=>(complex const&) const = default;
};
template <typename T> inline complex<T> conj(complex<T> const& x) { return x.conj(); }
template <typename T> inline T const real(complex<T> const& x) { return x.real(); }
template <typename T> inline T const imag(complex<T> const& x) { return x.imag(); }
template <typename T> inline T& real(complex<T>& x) { return x.real(); }
template <typename T> inline T& imag(complex<T>& x) { return x.imag(); }
template <typename T> inline constexpr complex<T> polar(T r, T theta) { return complex<T>::polar(r, theta); }

namespace math {
auto bpow(auto const& x, auto n, auto const& one, auto op) {
 if (n == 0) return one;
 auto t = bpow(x, n / 2, one, op);
 t = op(t, t);
 if (n % 2) t = op(t, x);
 return t;
}
auto bpow(auto x, auto n, auto ans) { return bpow(x, n, ans, std::multiplies{}); }
template <typename T> T bpow(T const& x, auto n) { return bpow(x, n, T(1)); }
inline constexpr auto inv2(auto x) {
 std::make_unsigned_t<decltype(x)> y = 1;
 while (y * x != 1) y*= 2 - x * y;
 return y;
}
} // namespace math
} // namespace cp_algo

namespace cp_algo::math::fft {
static constexpr size_t flen = 4;
using ftype = double;
using vftype = dx4;
using point = complex<ftype>;
using vpoint = complex<vftype>;
static constexpr vftype vz = {};
inline vpoint vi(vpoint const& r) { return {-imag(r), real(r)}; }

struct cvector {
 big_vector<vpoint> r;
 cvector(size_t n) {
  n = std::max(flen, std::bit_ceil(n));
  r.resize(n / flen);
  checkpoint("cvector create");
 }
 vpoint& at(size_t k) { return r[k / flen]; }
 vpoint at(size_t k) const { return r[k / flen]; }
 template <class pt = point> inline void set(size_t k, pt const& t) {
  if constexpr (std::is_same_v<pt, point>) {
   real(r[k / flen])[k % flen] = real(t);
   imag(r[k / flen])[k % flen] = imag(t);
  } else { at(k) = t; }
 }
 template <class pt = point> inline pt get(size_t k) const {
  if constexpr (std::is_same_v<pt, point>) {
   return {real(r[k / flen])[k % flen], imag(r[k / flen])[k % flen]};
  } else { return at(k); }
 }
 size_t size() const { return flen * r.size(); }
 static constexpr size_t eval_arg(size_t n) {
  if (n < pre_evals) return eval_args[n];
  return eval_arg(n / 2) | (n & 1) << (std::bit_width(n) - 1);
 }
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
 template <bool partial = true> void ifft() {
  size_t n = size();
  if constexpr (!partial) {
   point pi(0, 1);
   exec_on_evals<4>(n / 4, [&](size_t k, point rt) {
    k*= 4;
    point v1 = conj(rt);
    point v2 = v1 * v1;
    point v3 = v1 * v2;
    auto A = get(k);
    auto B = get(k + 1);
    auto C = get(k + 2);
    auto D = get(k + 3);
    set(k, (A + B) + (C + D));
    set(k + 2, ((A + B) - (C + D)) * v2);
    set(k + 1, ((A - B) - pi * (C - D)) * v1);
    set(k + 3, ((A - B) + pi * (C - D)) * v3);
   });
  }
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
  for (size_t k = 0; k < n; k+= flen) {
   if constexpr (partial) set(k, get<vpoint>(k)/= vz + ftype(n / flen));
   else set(k, get<vpoint>(k)/= vz + ftype(n));
  }
 }
 template <bool partial = true> void fft() {
  size_t n = size();
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
  if constexpr (!partial) {
   point pi(0, 1);
   exec_on_evals<4>(n / 4, [&](size_t k, point rt) {
    k*= 4;
    point v1 = rt;
    point v2 = v1 * v1;
    point v3 = v1 * v2;
    auto A = get(k);
    auto B = get(k + 1) * v1;
    auto C = get(k + 2) * v2;
    auto D = get(k + 3) * v3;
    set(k, (A + C) + (B + D));
    set(k + 1, (A + C) - (B + D));
    set(k + 2, (A - C) + pi * (B - D));
    set(k + 3, (A - C) - pi * (B - D));
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

template <bool simple = false> struct dft64 {
 big_vector<cvector> cv;
 static uint64_t factor, ifactor;
 static bool _init;
 static void init() {
  if (_init) return;
  _init = true;
  factor = simple ? 1 : random::rng();
  if (factor % 2 == 0) factor++;
  ifactor = inv2(factor);
 }
 static constexpr int parts = simple ? 1 : 4;
 dft64(auto const& a, size_t n): cv(parts, n) {
  init();
  uint64_t cur = 1, step = bpow(factor, n);
  for (size_t i = 0; i < std::min(std::size(a), n); i++) {
   auto split = [&](size_t i, uint64_t mul) -> std::array<int16_t, parts> {
    uint64_t x = i < std::size(a) ? a[i] * mul : 0;
    std::array<int16_t, parts> res;
    for (int z = 0; z < parts; z++) {
     res[z] = int16_t(x);
     x = (x >> 16) + (res[z] < 0);
    }
    return res;
   };
   auto re = split(i, cur);
   auto im = split(n + i, cur * step);
   for (int z = 0; z < parts; z++) {
    real(cv[z].at(i))[i % 4] = re[z];
    imag(cv[z].at(i))[i % 4] = im[z];
   }
   cur*= factor;
  }
  checkpoint("dft64 init");
  for (auto& x: cv) x.fft();
 }
 static void do_dot_iter(point rt, std::array<vpoint, parts>& B, std::array<vpoint, parts> const& A, std::array<vpoint, parts>& C) {
  for (size_t k = 0; k < parts; k++)
   for (size_t i = 0; i <= k; i++) C[k]+= A[i] * B[k - i];
  for (size_t k = 0; k < parts; k++) {
   real(B[k]) = rotate_right(real(B[k]));
   imag(B[k]) = rotate_right(imag(B[k]));
   auto bx = real(B[k])[0], by = imag(B[k])[0];
   real(B[k])[0] = bx * real(rt) - by * imag(rt);
   imag(B[k])[0] = bx * imag(rt) + by * real(rt);
  }
 }
 void dot(dft64 const& t) {
  size_t N = cv[0].size();
  cvector::exec_on_evals<1>(N / flen, [&](size_t k, point rt) {
   k*= flen;
   auto [A0x, A0y] = cv[0].at(k);
   vftype A1x{}, A1y{}, A2x{}, A2y{}, A3x{}, A3y{};
   if constexpr (!simple) {
    auto [_A1x, _A1y] = cv[1].at(k);
    auto [_A2x, _A2y] = cv[2].at(k);
    auto [_A3x, _A3y] = cv[3].at(k);
    A1x = _A1x; A1y = _A1y;
    A2x = _A2x; A2y = _A2y;
    A3x = _A3x; A3y = _A3y;
   }
   std::array<vpoint, parts> A, B, C;
   if constexpr (simple) {
    B = {t.cv[0].at(k)};
    C = {vz};
   } else {
    B = {t.cv[0].at(k), t.cv[1].at(k), t.cv[2].at(k), t.cv[3].at(k)};
    C = {vz, vz, vz, vz};
   }
   for (size_t i = 0; i < flen; i++) {
    if constexpr (simple) {
     A = {vpoint{vz + A0x[i], vz + A0y[i]}};
    } else {
     A = {vpoint{vz + A0x[i], vz + A0y[i]}, vpoint{vz + A1x[i], vz + A1y[i]}, vpoint{vz + A2x[i], vz + A2y[i]}, vpoint{vz + A3x[i], vz + A3y[i]}};
    }
    do_dot_iter(rt, B, A, C);
   }
   cv[0].at(k) = C[0];
   if constexpr (!simple) {
    cv[1].at(k) = C[1];
    cv[2].at(k) = C[2];
    cv[3].at(k) = C[3];
   }
  });
  checkpoint("dot");
  for (auto& x: cv) x.ifft();
 }
 void recover_mod(auto& res, size_t k) {
  size_t n = cv[0].size();
  uint64_t cur = 1, step = bpow(ifactor, n);
  for (size_t i = 0; i < std::min(k, n); i++) {
   std::array<ftype, parts> re, im;
   if constexpr (simple) {
    re = {real(cv[0].get(i))};
    im = {imag(cv[0].get(i))};
   } else {
    re = {real(cv[0].get(i)), real(cv[1].get(i)), real(cv[2].get(i)), real(cv[3].get(i))};
    im = {imag(cv[0].get(i)), imag(cv[1].get(i)), imag(cv[2].get(i)), imag(cv[3].get(i))};
   }
   auto set_i = [&](size_t i, auto& x, auto mul) {
    if (i >= k) return;
    if constexpr (simple) {
     res[i] = (uint64_t)(int64_t) llround(x[0]);
    } else {
     res[i] = (uint64_t)(int64_t) llround(x[0]) + ((uint64_t)(int64_t) llround(x[1]) << 16) + ((uint64_t)(int64_t) llround(x[2]) << 32) + ((uint64_t)(int64_t) llround(x[3]) << 48);
    }
    res[i]*= mul;
   };
   set_i(i, re, cur);
   set_i(n + i, im, cur * step);
   cur*= ifactor;
  }
  cp_algo::checkpoint("recover mod");
 }
};
template <bool simple> uint64_t dft64<simple>::factor = 1;
template <bool simple> uint64_t dft64<simple>::ifactor = 1;
template <bool simple> bool dft64<simple>::_init = false;

template <bool simple = false> void conv64(auto& a, auto const& b) {
 size_t n = a.size(), m = b.size();
 size_t N = std::max(flen, std::bit_ceil(n + m - 1) / 2);
 dft64<simple> A(a, N), B(b, N);
 A.dot(B);
 a.resize(n + m - 1);
 A.recover_mod(a, n + m - 1);
}
} // namespace cp_algo::math::fft

struct Conv {
 static vector<u64> run(const vector<u64>& a, const vector<u64>& b) {
  if (a.empty() || b.empty()) return {};
  std::vector<uint64_t, cp_algo::big_alloc<uint64_t>> aa(a.begin(), a.end()), bb(b.begin(), b.end());
  cp_algo::math::fft::conv64(aa, bb);
  return vector<u64>(aa.begin(), aa.end());
 }
};
