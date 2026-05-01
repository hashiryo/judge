#pragma once
// =============================================================================
// Source: yosupo "Characteristic Polynomial" 提出 285790 (cp-algo).
//   https://judge.yosupo.jp/submission/285790
//   方式: cp-algo の linalg::frobenius_form(A) で Frobenius 正規形を求め、
//   その不変因子 (各ブロックの最小多項式) を掛け合わせて A の特性多項式とする。
//   内部は modint + linalg::matrix + math::poly (FFT, Euclid GCD, half_gcd 等)
//   の重い template 連鎖。FFT は AVX2 (_mm256_mul_epu32) montgomery + 自前
//   complex<double> + cvector で実装され、N≤500 でも N^ω 系オーダで動く。
// 抽出方針:
//   - blazingio I/O は使わない (base.cpp が cin で済ませる)
//   - <experimental/simd> include を削除 (libstdc++ stdfloat::float16_t と
//     arm_neon.h ::float16_t の衝突回避; cp-algo 自身は stdx::* を使っていない)
//   - mmap path を Linux 限定に (MADV_HUGEPAGE/POPULATE_WRITE が macOS にない)
//   - __builtin_shuffle を clang では __builtin_shufflevector に分岐
//   - #line マーカー除去
//   - CharPoly::run wrapper で vector<vector<u32>> → matrix<base> に詰めて
//     frobenius_form を呼び、各ブロックの polyn を掛け合わせる
// ライセンス: cp-algo (元実装に従う)
// =============================================================================

#pragma GCC optimize("Ofast,unroll-loops")
#define CP_ALGO_MAXN (1 << 10)
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2,bmi,bmi2")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#elif defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "_common.hpp"

#include <ranges>
#include <numeric>









#include <optional>

    namespace cp_algo::math
{
 // a * x + b
 template <typename base>
 struct lin
 {
  base a = 1, b = 0;
  std::optional<base> c;
  lin() {}
  lin(base b) : a(0), b(b) {}
  lin(base a, base b) : a(a), b(b) {}
  lin(base a, base b, base _c) : a(a), b(b), c(_c) {}

  // polynomial product modulo x^2 - c
  lin operator*(const lin &t)
  {
   assert(c && t.c && *c == *t.c);
   return {a * t.b + b * t.a, b * t.b + a * t.a * (*c), *c};
  }

  // a * (t.a * x + t.b) + b
  lin apply(lin const &t) const
  {
   return {a * t.a, a * t.b + b};
  }

  void prepend(lin const &t)
  {
   *this = t.apply(*this);
  }

  base eval(base x) const
  {
   return a * x + b;
  }
 };

 // (ax+b) / (cx+d)
 template <typename base>
 struct linfrac
 {
  base a, b, c, d;
  linfrac() : a(1), b(0), c(0), d(1) {}       // x, identity for composition
  linfrac(base a) : a(a), b(1), c(1), d(0) {} // a + 1/x, for continued fractions
  linfrac(base a, base b, base c, base d) : a(a), b(b), c(c), d(d) {}

  // composition of two linfracs
  linfrac operator*(linfrac t) const
  {
   return t.prepend(linfrac(*this));
  }

  linfrac operator-() const
  {
   return {-a, -b, -c, -d};
  }

  linfrac adj() const
  {
   return {d, -b, -c, a};
  }

  linfrac &prepend(linfrac const &t)
  {
   t.apply(a, c);
   t.apply(b, d);
   return *this;
  }

  // apply linfrac to A/B
  void apply(base &A, base &B) const
  {
   std::tie(A, B) = std::pair{a * A + b * B, c * A + d * B};
  }
 };
}








namespace cp_algo::math
{
#ifdef CP_ALGO_MAXN
 const int maxn = CP_ALGO_MAXN;
#else
 const int maxn = 1 << 19;
#endif
 const int magic = 64; // threshold for sizes to run the naive algo

 auto bpow(auto const &x, auto n, auto const &one, auto op)
 {
  if (n == 0)
  {
   return one;
  }
  else
  {
   auto t = bpow(x, n / 2, one, op);
   t = op(t, t);
   if (n % 2)
   {
    t = op(t, x);
   }
   return t;
  }
 }
 auto bpow(auto x, auto n, auto ans)
 {
  return bpow(x, n, ans, std::multiplies{});
 }
 template <typename T>
 T bpow(T const &x, auto n)
 {
  return bpow(x, n, T(1));
 }
}


namespace cp_algo::math
{

 template <typename modint, typename _Int>
 struct modint_base
 {
  using Int = _Int;
  using UInt = std::make_unsigned_t<Int>;
  static constexpr size_t bits = sizeof(Int) * 8;
  using Int2 = std::conditional_t<bits <= 32, int64_t, __int128_t>;
  using UInt2 = std::conditional_t<bits <= 32, uint64_t, __uint128_t>;
  static Int mod()
  {
   return modint::mod();
  }
  static Int remod()
  {
   return modint::remod();
  }
  static UInt2 modmod()
  {
   return UInt2(mod()) * mod();
  }
  modint_base() = default;
  modint_base(Int2 rr)
  {
   to_modint().setr(UInt((rr + modmod()) % mod()));
  }
  modint inv() const
  {
   return bpow(to_modint(), mod() - 2);
  }
  modint operator-() const
  {
   modint neg;
   neg.r = std::min(-r, remod() - r);
   return neg;
  }
  modint &operator/=(const modint &t)
  {
   return to_modint() *= t.inv();
  }
  modint &operator*=(const modint &t)
  {
   r = UInt(UInt2(r) * t.r % mod());
   return to_modint();
  }
  modint &operator+=(const modint &t)
  {
   r += t.r;
   r = std::min(r, r - remod());
   return to_modint();
  }
  modint &operator-=(const modint &t)
  {
   r -= t.r;
   r = std::min(r, r + remod());
   return to_modint();
  }
  modint operator+(const modint &t) const { return modint(to_modint()) += t; }
  modint operator-(const modint &t) const { return modint(to_modint()) -= t; }
  modint operator*(const modint &t) const { return modint(to_modint()) *= t; }
  modint operator/(const modint &t) const { return modint(to_modint()) /= t; }
  // Why <=> doesn't work?..
  auto operator==(const modint &t) const { return to_modint().getr() == t.getr(); }
  auto operator!=(const modint &t) const { return to_modint().getr() != t.getr(); }
  auto operator<=(const modint &t) const { return to_modint().getr() <= t.getr(); }
  auto operator>=(const modint &t) const { return to_modint().getr() >= t.getr(); }
  auto operator<(const modint &t) const { return to_modint().getr() < t.getr(); }
  auto operator>(const modint &t) const { return to_modint().getr() > t.getr(); }
  Int rem() const
  {
   UInt R = to_modint().getr();
   return R - (R > (UInt)mod() / 2) * mod();
  }
  void setr(UInt rr)
  {
   r = rr;
  }
  UInt getr() const
  {
   return r;
  }

  // Only use these if you really know what you're doing!
  static UInt modmod8() { return UInt(8 * modmod()); }
  void add_unsafe(UInt t) { r += t; }
  void pseudonormalize() { r = std::min(r, r - modmod8()); }
  modint const &normalize()
  {
   if (r >= (UInt)mod())
   {
    r %= mod();
   }
   return to_modint();
  }
  void setr_direct(UInt rr) { r = rr; }
  UInt getr_direct() const { return r; }

 protected:
  UInt r;

 private:
  modint &to_modint() { return static_cast<modint &>(*this); }
  modint const &to_modint() const { return static_cast<modint const &>(*this); }
 };
 template <typename modint>
 concept modint_type = std::is_base_of_v<modint_base<modint, typename modint::Int>, modint>;
 template <modint_type modint>
 decltype(std::cin) &operator>>(decltype(std::cin) &in, modint &x)
 {
  typename modint::UInt r;
  auto &res = in >> r;
  x.setr(r);
  return res;
 }
 template <modint_type modint>
 decltype(std::cout) &operator<<(decltype(std::cout) &out, modint const &x)
 {
  return out << x.getr();
 }

 template <auto m>
 struct modint : modint_base<modint<m>, decltype(m)>
 {
  using Base = modint_base<modint<m>, decltype(m)>;
  using Base::Base;
  static constexpr Base::Int mod() { return m; }
  static constexpr Base::UInt remod() { return m; }
  auto getr() const { return Base::r; }
 };

 inline constexpr auto inv2(auto x)
 {
  assert(x % 2);
  std::make_unsigned_t<decltype(x)> y = 1;
  while (y * x != 1)
  {
   y *= 2 - x * y;
  }
  return y;
 }

 template <typename Int = int64_t>
 struct dynamic_modint : modint_base<dynamic_modint<Int>, Int>
 {
  using Base = modint_base<dynamic_modint<Int>, Int>;
  using Base::Base;

  static Base::UInt m_reduce(Base::UInt2 ab)
  {
   if (mod() % 2 == 0) [[unlikely]]
   {
    return typename Base::UInt(ab % mod());
   }
   else
   {
    typename Base::UInt2 m = typename Base::UInt(ab) * imod();
    return typename Base::UInt((ab + m * mod()) >> Base::bits);
   }
  }
  static Base::UInt m_transform(Base::UInt a)
  {
   if (mod() % 2 == 0) [[unlikely]]
   {
    return a;
   }
   else
   {
    return m_reduce(a * pw128());
   }
  }
  dynamic_modint &operator*=(const dynamic_modint &t)
  {
   Base::r = m_reduce(typename Base::UInt2(Base::r) * t.r);
   return *this;
  }
  void setr(Base::UInt rr)
  {
   Base::r = m_transform(rr);
  }
  Base::UInt getr() const
  {
   typename Base::UInt res = m_reduce(Base::r);
   return std::min(res, res - mod());
  }
  static Int mod() { return m; }
  static Int remod() { return 2 * m; }
  static Base::UInt imod() { return im; }
  static Base::UInt2 pw128() { return r2; }
  static void switch_mod(Int nm)
  {
   m = nm;
   im = m % 2 ? inv2(-m) : 0;
   r2 = static_cast<Base::UInt>(static_cast<Base::UInt2>(-1) % m + 1);
  }

  // Wrapper for temp switching
  auto static with_mod(Int tmp, auto callback)
  {
   struct scoped
   {
    Int prev = mod();
    ~scoped() { switch_mod(prev); }
   } _;
   switch_mod(tmp);
   return callback();
  }

 private:
  static thread_local Int m;
  static thread_local Base::UInt im, r2;
 };
 template <typename Int>
 Int thread_local dynamic_modint<Int>::m = 1;
 template <typename Int>
 dynamic_modint<Int>::Base::UInt thread_local dynamic_modint<Int>::im = -1;
 template <typename Int>
 dynamic_modint<Int>::Base::UInt thread_local dynamic_modint<Int>::r2 = 0;
}




namespace cp_algo
{
 std::map<std::string, double> checkpoints;
 template <bool final = false>
 void checkpoint([[maybe_unused]] std::string const &msg = "")
 {
#ifdef CP_ALGO_CHECKPOINT
  static double last = 0;
  double now = (double)clock() / CLOCKS_PER_SEC;
  double delta = now - last;
  last = now;
  if (msg.size() && !final)
  {
   checkpoints[msg] += delta;
  }
  if (final)
  {
   for (auto const &[key, value] : checkpoints)
   {
    std::cerr << key << ": " << value * 1000 << " ms\n";
   }
   std::cerr << "Total: " << now * 1000 << " ms\n";
  }
#endif
 }
}




namespace cp_algo::random
{
 uint64_t rng()
 {
  static std::mt19937_64 rng(
      std::chrono::steady_clock::now().time_since_epoch().count());
  return rng();
 }
}







namespace cp_algo
{
 template <typename T, size_t len>
 using simd [[gnu::vector_size(len * sizeof(T))]] = T;
 using i64x4 = simd<int64_t, 4>;
 using u64x4 = simd<uint64_t, 4>;
 using u32x8 = simd<uint32_t, 8>;
 using i32x4 = simd<int32_t, 4>;
 using u32x4 = simd<uint32_t, 4>;
 using dx4 = simd<double, 4>;

 [[gnu::always_inline]] inline dx4 abs(dx4 a)
 {
  return a < 0 ? -a : a;
 }

 // https://stackoverflow.com/a/77376595
 // works for ints in (-2^51, 2^51)
 static constexpr dx4 magic = dx4() + (3ULL << 51);
 [[gnu::always_inline]] inline i64x4 lround(dx4 x)
 {
  return i64x4(x + magic) - i64x4(magic);
 }
 [[gnu::always_inline]] inline dx4 to_double(i64x4 x)
 {
  return dx4(x + i64x4(magic)) - magic;
 }

 [[gnu::always_inline]] inline dx4 round(dx4 a)
 {
  return dx4{
      std::nearbyint(a[0]),
      std::nearbyint(a[1]),
      std::nearbyint(a[2]),
      std::nearbyint(a[3])};
 }

 [[gnu::always_inline]] inline u64x4 montgomery_reduce(u64x4 x, u64x4 mod, u64x4 imod)
 {
  auto x_ninv = u64x4(u32x8(x) * u32x8(imod));
#ifdef __AVX2__
  x += u64x4(_mm256_mul_epu32(__m256i(x_ninv), __m256i(mod)));
#else
  x += x_ninv * mod;
#endif
  return x >> 32;
 }

 [[gnu::always_inline]] inline u64x4 montgomery_mul(u64x4 x, u64x4 y, u64x4 mod, u64x4 imod)
 {
#ifdef __AVX2__
  return montgomery_reduce(u64x4(_mm256_mul_epu32(__m256i(x), __m256i(y))), mod, imod);
#else
  return montgomery_reduce(x * y, mod, imod);
#endif
 }

 [[gnu::always_inline]] inline dx4 rotate_right(dx4 x)
 {
#if defined(__clang__)
  return __builtin_shufflevector(x, x, 3, 0, 1, 2);
#else
  static constexpr u64x4 shuffler = {3, 0, 1, 2};
  return __builtin_shuffle(x, shuffler);
#endif
 }

 template <std::size_t Align = 32>
 [[gnu::always_inline]] inline bool is_aligned(const auto *p) noexcept
 {
  return (reinterpret_cast<std::uintptr_t>(p) % Align) == 0;
 }

 template <class Target>
 [[gnu::always_inline]] inline Target &vector_cast(auto &&p)
 {
  return *reinterpret_cast<Target *>(std::assume_aligned<alignof(Target)>(&p));
 }
}




namespace cp_algo
{
 // Custom implementation, since std::complex is UB on non-floating types
 template <typename T>
 struct complex
 {
  using value_type = T;
  T x, y;
  constexpr complex() : x(), y() {}
  constexpr complex(T x) : x(x), y() {}
  constexpr complex(T x, T y) : x(x), y(y) {}
  complex &operator*=(T t)
  {
   x *= t;
   y *= t;
   return *this;
  }
  complex &operator/=(T t)
  {
   x /= t;
   y /= t;
   return *this;
  }
  complex operator*(T t) const { return complex(*this) *= t; }
  complex operator/(T t) const { return complex(*this) /= t; }
  complex &operator+=(complex t)
  {
   x += t.x;
   y += t.y;
   return *this;
  }
  complex &operator-=(complex t)
  {
   x -= t.x;
   y -= t.y;
   return *this;
  }
  complex operator*(complex t) const { return {x * t.x - y * t.y, x * t.y + y * t.x}; }
  complex operator/(complex t) const { return *this * t.conj() / t.norm(); }
  complex operator+(complex t) const { return complex(*this) += t; }
  complex operator-(complex t) const { return complex(*this) -= t; }
  complex &operator*=(complex t) { return *this = *this * t; }
  complex &operator/=(complex t) { return *this = *this / t; }
  complex operator-() const { return {-x, -y}; }
  complex conj() const { return {x, -y}; }
  T norm() const { return x * x + y * y; }
  T abs() const { return std::sqrt(norm()); }
  T const real() const { return x; }
  T const imag() const { return y; }
  T &real() { return x; }
  T &imag() { return y; }
  static constexpr complex polar(T r, T theta) { return {T(r * cos(theta)), T(r * sin(theta))}; }
  auto operator<=>(complex const &t) const = default;
 };
 template <typename T>
 complex<T> operator*(auto x, complex<T> y) { return y *= x; }
 template <typename T>
 complex<T> conj(complex<T> x) { return x.conj(); }
 template <typename T>
 T norm(complex<T> x) { return x.norm(); }
 template <typename T>
 T abs(complex<T> x) { return x.abs(); }
 template <typename T>
 T &real(complex<T> &x) { return x.real(); }
 template <typename T>
 T &imag(complex<T> &x) { return x.imag(); }
 template <typename T>
 T const real(complex<T> const &x) { return x.real(); }
 template <typename T>
 T const imag(complex<T> const &x) { return x.imag(); }
 template <typename T>
 constexpr complex<T> polar(T r, T theta)
 {
  return complex<T>::polar(r, theta);
 }
 template <typename T>
 std::ostream &operator<<(std::ostream &out, complex<T> x)
 {
  return out << x.real() << ' ' << x.imag();
 }
}





// macOS は MADV_HUGEPAGE / MADV_POPULATE_WRITE が無いので Linux 限定にする。
#if defined(__linux__)
#define CP_ALGO_USE_MMAP 1
#include <sys/mman.h>
#else
#define CP_ALGO_USE_MMAP 0
#endif

namespace cp_algo
{
 template <typename T, std::size_t Align = 32>
 class big_alloc
 {
  static_assert(Align >= alignof(void *), "Align must be at least pointer-size");
  static_assert(std::popcount(Align) == 1, "Align must be a power of two");

 public:
  using value_type = T;
  template <class U>
  struct rebind
  {
   using other = big_alloc<U, Align>;
  };
  constexpr bool operator==(const big_alloc &) const = default;
  constexpr bool operator!=(const big_alloc &) const = default;

  big_alloc() noexcept = default;
  template <typename U, std::size_t A>
  big_alloc(const big_alloc<U, A> &) noexcept {}

  [[nodiscard]] T *allocate(std::size_t n)
  {
   std::size_t padded = round_up(n * sizeof(T));
   std::size_t align = std::max<std::size_t>(alignof(T), Align);
#if CP_ALGO_USE_MMAP
   if (padded >= MEGABYTE)
   {
    void *raw = mmap(nullptr, padded,
                     PROT_READ | PROT_WRITE,
                     MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    madvise(raw, padded, MADV_HUGEPAGE);
    madvise(raw, padded, MADV_POPULATE_WRITE);
    return static_cast<T *>(raw);
   }
#endif
   return static_cast<T *>(::operator new(padded, std::align_val_t(align)));
  }

  void deallocate(T *p, std::size_t n) noexcept
  {
   if (!p)
    return;
   std::size_t padded = round_up(n * sizeof(T));
   std::size_t align = std::max<std::size_t>(alignof(T), Align);
#if CP_ALGO_USE_MMAP
   if (padded >= MEGABYTE)
   {
    munmap(p, padded);
    return;
   }
#endif
   ::operator delete(p, padded, std::align_val_t(align));
  }

 private:
  static constexpr std::size_t MEGABYTE = 1 << 20;
  static constexpr std::size_t round_up(std::size_t x) noexcept
  {
   return (x + Align - 1) / Align * Align;
  }
 };
}


#include <ranges>
#include <bit>

// std::experimental::simd は include を消したので alias も削除 (cp-algo は使っていない)
namespace cp_algo::math::fft
{
 static constexpr size_t flen = 4;
 using ftype = double;
 using vftype = dx4;
 using point = complex<ftype>;
 using vpoint = complex<vftype>;
 static constexpr vftype vz = {};
 vpoint vi(vpoint const &r)
 {
  return {-imag(r), real(r)};
 }

 struct cvector
 {
  std::vector<vpoint, big_alloc<vpoint>> r;
  cvector(size_t n)
  {
   n = std::max(flen, std::bit_ceil(n));
   r.resize(n / flen);
   checkpoint("cvector create");
  }

  vpoint &at(size_t k) { return r[k / flen]; }
  vpoint at(size_t k) const { return r[k / flen]; }
  template <class pt = point>
  void set(size_t k, pt t)
  {
   if constexpr (std::is_same_v<pt, point>)
   {
    real(r[k / flen])[k % flen] = real(t);
    imag(r[k / flen])[k % flen] = imag(t);
   }
   else
   {
    at(k) = t;
   }
  }
  template <class pt = point>
  pt get(size_t k) const
  {
   if constexpr (std::is_same_v<pt, point>)
   {
    return {real(r[k / flen])[k % flen], imag(r[k / flen])[k % flen]};
   }
   else
   {
    return at(k);
   }
  }

  size_t size() const
  {
   return flen * r.size();
  }
  static constexpr size_t eval_arg(size_t n)
  {
   if (n < pre_evals)
   {
    return eval_args[n];
   }
   else
   {
    return eval_arg(n / 2) | (n & 1) << (std::bit_width(n) - 1);
   }
  }
  static constexpr point eval_point(size_t n)
  {
   if (n % 2)
   {
    return -eval_point(n - 1);
   }
   else if (n % 4)
   {
    return eval_point(n - 2) * point(0, 1);
   }
   else if (n / 4 < pre_evals)
   {
    return evalp[n / 4];
   }
   else
   {
    return polar<ftype>(1., std::numbers::pi / (ftype)std::bit_floor(n) * (ftype)eval_arg(n));
   }
  }
  // libc++ では cos/sin が非 constexpr なので const inline に格下げ。
  static inline const std::array<point, 32> roots = []()
  {
   std::array<point, 32> res;
   for (size_t i = 2; i < 32; i++)
   {
    res[i] = polar<ftype>(1., std::numbers::pi / (1ull << (i - 2)));
   }
   return res;
  }();
  static point root(size_t n)
  {
   return roots[std::bit_width(n)];
  }
  template <int step>
  static void exec_on_eval(size_t n, size_t k, auto &&callback)
  {
   callback(k, root(4 * step * n) * eval_point(step * k));
  }
  template <int step>
  static void exec_on_evals(size_t n, auto &&callback)
  {
   point factor = root(4 * step * n);
   for (size_t i = 0; i < n; i++)
   {
    callback(i, factor * eval_point(step * i));
   }
  }

  void dot(cvector const &t)
  {
   size_t n = this->size();
   exec_on_evals<1>(n / flen, [&](size_t k, point rt)
                    {
                k *= flen;
                auto [Ax, Ay] = at(k);
                auto Bv = t.at(k);
                vpoint res = vz;
                for (size_t i = 0; i < flen; i++) {
                    res += vpoint(vz + Ax[i], vz + Ay[i]) * Bv;
                    real(Bv) = rotate_right(real(Bv));
                    imag(Bv) = rotate_right(imag(Bv));
                    auto x = real(Bv)[0], y = imag(Bv)[0];
                    real(Bv)[0] = x * real(rt) - y * imag(rt);
                    imag(Bv)[0] = x * imag(rt) + y * real(rt);
                }
                set(k, res); });
   checkpoint("dot");
  }

  void ifft()
  {
   size_t n = size();
   bool parity = std::countr_zero(n) % 2;
   if (parity)
   {
    exec_on_evals<2>(n / (2 * flen), [&](size_t k, point rt)
                     {
                    k *= 2 * flen;
                    vpoint cvrt = {vz + real(rt), vz - imag(rt)};
                    auto B = at(k) - at(k + flen);
                    at(k) += at(k + flen);
                    at(k + flen) = B * cvrt; });
   }

   for (size_t leaf = 3 * flen; leaf < n; leaf += 4 * flen)
   {
    size_t level = std::countr_one(leaf + 3);
    for (size_t lvl = 4 + parity; lvl <= level; lvl += 2)
    {
     size_t i = (1 << lvl) / 4;
     exec_on_eval<4>(n >> lvl, leaf >> lvl, [&](size_t k, point rt)
                     {
                        k <<= lvl;
                        vpoint v1 = {vz + real(rt), vz - imag(rt)};
                        vpoint v2 = v1 * v1;
                        vpoint v3 = v1 * v2;
                        for(size_t j = k; j < k + i; j += flen) {
                            auto A = at(j);
                            auto B = at(j + i);
                            auto C = at(j + 2 * i);
                            auto D = at(j + 3 * i);
                            at(j) = ((A + B) + (C + D));
                            at(j + 2 * i) = ((A + B) - (C + D)) * v2;
                            at(j +     i) = ((A - B) - vi(C - D)) * v1;
                            at(j + 3 * i) = ((A - B) + vi(C - D)) * v3;
                        } });
    }
   }
   checkpoint("ifft");
   for (size_t k = 0; k < n; k += flen)
   {
    set(k, get<vpoint>(k) /= vz + (ftype)(n / flen));
   }
  }
  void fft()
  {
   size_t n = size();
   bool parity = std::countr_zero(n) % 2;
   for (size_t leaf = 0; leaf < n; leaf += 4 * flen)
   {
    size_t level = std::countr_zero(n + leaf);
    level -= level % 2 != parity;
    for (size_t lvl = level; lvl >= 4; lvl -= 2)
    {
     size_t i = (1 << lvl) / 4;
     exec_on_eval<4>(n >> lvl, leaf >> lvl, [&](size_t k, point rt)
                     {
                        k <<= lvl;
                        vpoint v1 = {vz + real(rt), vz + imag(rt)};
                        vpoint v2 = v1 * v1;
                        vpoint v3 = v1 * v2;
                        for(size_t j = k; j < k + i; j += flen) {
                            auto A = at(j);
                            auto B = at(j + i) * v1;
                            auto C = at(j + 2 * i) * v2;
                            auto D = at(j + 3 * i) * v3;
                            at(j)         = (A + C) + (B + D);
                            at(j + i)     = (A + C) - (B + D);
                            at(j + 2 * i) = (A - C) + vi(B - D);
                            at(j + 3 * i) = (A - C) - vi(B - D);
                        } });
    }
   }
   if (parity)
   {
    exec_on_evals<2>(n / (2 * flen), [&](size_t k, point rt)
                     {
                    k *= 2 * flen;
                    vpoint vrt = {vz + real(rt), vz + imag(rt)};
                    auto t = at(k + flen) * vrt;
                    at(k + flen) = at(k) - t;
                    at(k) += t; });
   }
   checkpoint("fft");
  }
  static constexpr size_t pre_evals = 1 << 16;
  static const std::array<size_t, pre_evals> eval_args;
  static const std::array<point, pre_evals> evalp;
 };

 const std::array<size_t, cvector::pre_evals> cvector::eval_args = []()
 {
  std::array<size_t, pre_evals> res = {};
  for (size_t i = 1; i < pre_evals; i++)
  {
   res[i] = res[i >> 1] | (i & 1) << (std::bit_width(i) - 1);
  }
  return res;
 }();
 const std::array<point, cvector::pre_evals> cvector::evalp = []()
 {
  std::array<point, pre_evals> res = {};
  res[0] = 1;
  for (size_t n = 1; n < pre_evals; n++)
  {
   res[n] = polar<ftype>(1., std::numbers::pi * ftype(eval_args[n]) / ftype(4 * std::bit_floor(n)));
  }
  return res;
 }();
}


namespace cp_algo::math::fft
{
 template <modint_type base>
 struct dft
 {
  cvector A, B;
  static base factor, ifactor;
  using Int2 = base::Int2;
  static bool _init;
  static int split()
  {
   static const int splt = int(std::sqrt(base::mod())) + 1;
   return splt;
  }
  static u64x4 mod, imod;

  static void init()
  {
   if (!_init)
   {
    factor = 1 + random::rng() % (base::mod() - 1);
    ifactor = base(1) / factor;
    mod = u64x4() + base::mod();
    imod = u64x4() + inv2(-base::mod());
    _init = true;
   }
  }

  dft(auto const &a, size_t n) : A(n), B(n)
  {
   init();
   base b2x32 = bpow(base(2), 32);
   u64x4 cur = {
       (bpow(factor, 1) * b2x32).getr(),
       (bpow(factor, 2) * b2x32).getr(),
       (bpow(factor, 3) * b2x32).getr(),
       (bpow(factor, 4) * b2x32).getr()};
   u64x4 step4 = u64x4{} + (bpow(factor, 4) * b2x32).getr();
   u64x4 stepn = u64x4{} + (bpow(factor, n) * b2x32).getr();
   for (size_t i = 0; i < std::min(n, std::size(a)); i += flen)
   {
    auto splt = [&](size_t i, auto mul)
    {
     if (i >= std::size(a))
     {
      return std::pair{vftype(), vftype()};
     }
     u64x4 au = {
         i < std::size(a) ? a[i].getr() : 0,
         i + 1 < std::size(a) ? a[i + 1].getr() : 0,
         i + 2 < std::size(a) ? a[i + 2].getr() : 0,
         i + 3 < std::size(a) ? a[i + 3].getr() : 0};
     au = montgomery_mul(au, mul, mod, imod);
     au = au >= base::mod() ? au - base::mod() : au;
     auto ai = to_double(i64x4(au >= base::mod() / 2 ? au - base::mod() : au));
     auto quo = round(ai / split());
     return std::pair{ai - quo * split(), quo};
    };
    auto [rai, qai] = splt(i, cur);
    auto [rani, qani] = splt(n + i, montgomery_mul(cur, stepn, mod, imod));
    A.at(i) = vpoint(rai, rani);
    B.at(i) = vpoint(qai, qani);
    cur = montgomery_mul(cur, step4, mod, imod);
   }
   checkpoint("dft init");
   if (n)
   {
    A.fft();
    B.fft();
   }
  }

  void dot(auto &&C, auto const &D)
  {
   cvector::exec_on_evals<1>(A.size() / flen, [&](size_t k, point rt)
                             {
                k *= flen;
                auto [Ax, Ay] = A.at(k);
                auto [Bx, By] = B.at(k);
                vpoint AC, AD, BC, BD;
                AC = AD = BC = BD = vz;
                auto Cv = C.at(k), Dv = D.at(k);
                for (size_t i = 0; i < flen; i++) {
                    vpoint Av = {vz + Ax[i], vz + Ay[i]}, Bv = {vz + Bx[i], vz + By[i]};
                    AC += Av * Cv; AD += Av * Dv;
                    BC += Bv * Cv; BD += Bv * Dv;
                    real(Cv) = rotate_right(real(Cv));
                    imag(Cv) = rotate_right(imag(Cv));
                    real(Dv) = rotate_right(real(Dv));
                    imag(Dv) = rotate_right(imag(Dv));
                    auto cx = real(Cv)[0], cy = imag(Cv)[0];
                    auto dx = real(Dv)[0], dy = imag(Dv)[0];
                    real(Cv)[0] = cx * real(rt) - cy * imag(rt);
                    imag(Cv)[0] = cx * imag(rt) + cy * real(rt);
                    real(Dv)[0] = dx * real(rt) - dy * imag(rt);
                    imag(Dv)[0] = dx * imag(rt) + dy * real(rt);
                }
                A.at(k) = AC;
                C.at(k) = AD + BC;
                B.at(k) = BD; });
   checkpoint("dot");
  }

  void recover_mod(auto &&C, auto &res, size_t k)
  {
   size_t check = (k + flen - 1) / flen * flen;
   assert(res.size() >= check);
   size_t n = A.size();
   auto const splitsplit = base(split() * split()).getr();
   base b2x32 = bpow(base(2), 32);
   base b2x64 = bpow(base(2), 64);
   u64x4 cur = {
       (bpow(ifactor, 2) * b2x64).getr(),
       (bpow(ifactor, 3) * b2x64).getr(),
       (bpow(ifactor, 4) * b2x64).getr(),
       (bpow(ifactor, 5) * b2x64).getr()};
   u64x4 step4 = u64x4{} + (bpow(ifactor, 4) * b2x32).getr();
   u64x4 stepn = u64x4{} + (bpow(ifactor, n) * b2x32).getr();
   for (size_t i = 0; i < std::min(n, k); i += flen)
   {
    auto [Ax, Ay] = A.at(i);
    auto [Bx, By] = B.at(i);
    auto [Cx, Cy] = C.at(i);
    auto set_i = [&](size_t i, auto A, auto B, auto C, auto mul)
    {
     auto A0 = lround(A), A1 = lround(C), A2 = lround(B);
     // clang は i64x4 × uint64_t の scalar→vector 拡張を許さない (truncation
     // 警告)。明示的に scalar → i64x4 に揃える。
     auto Ai = A0 + A1 * (int64_t) split() + A2 * (int64_t) splitsplit + (i64x4) (u64x4{} + (uint64_t) base::modmod());
     auto Au = montgomery_reduce(u64x4(Ai), mod, imod);
     Au = montgomery_mul(Au, mul, mod, imod);
     Au = Au >= base::mod() ? Au - base::mod() : Au;
     for (size_t j = 0; j < flen; j++)
     {
      res[i + j].setr(typename base::UInt(Au[j]));
     }
    };
    set_i(i, Ax, Bx, Cx, cur);
    if (i + n < k)
    {
     set_i(i + n, Ay, By, Cy, montgomery_mul(cur, stepn, mod, imod));
    }
    cur = montgomery_mul(cur, step4, mod, imod);
   }
   checkpoint("recover mod");
  }

  void mul(auto &&C, auto const &D, auto &res, size_t k)
  {
   assert(A.size() == C.size());
   size_t n = A.size();
   if (!n)
   {
    res = {};
    return;
   }
   dot(C, D);
   A.ifft();
   B.ifft();
   C.ifft();
   recover_mod(C, res, k);
  }
  void mul_inplace(auto &&B, auto &res, size_t k)
  {
   mul(B.A, B.B, res, k);
  }
  void mul(auto const &B, auto &res, size_t k)
  {
   mul(cvector(B.A), B.B, res, k);
  }
  std::vector<base, big_alloc<base>> operator*=(dft &B)
  {
   std::vector<base, big_alloc<base>> res(2 * A.size());
   mul_inplace(B, res, 2 * A.size());
   return res;
  }
  std::vector<base, big_alloc<base>> operator*=(dft const &B)
  {
   std::vector<base, big_alloc<base>> res(2 * A.size());
   mul(B, res, 2 * A.size());
   return res;
  }
  auto operator*(dft const &B) const
  {
   return dft(*this) *= B;
  }

  point operator[](int i) const { return A.get(i); }
 };
 template <modint_type base>
 base dft<base>::factor = 1;
 template <modint_type base>
 base dft<base>::ifactor = 1;
 template <modint_type base>
 bool dft<base>::_init = false;
 template <modint_type base>
 u64x4 dft<base>::mod = {};
 template <modint_type base>
 u64x4 dft<base>::imod = {};

 void mul_slow(auto &a, auto const &b, size_t k)
 {
  if (std::empty(a) || std::empty(b))
  {
   a.clear();
  }
  else
  {
   size_t n = std::min(k, std::size(a));
   size_t m = std::min(k, std::size(b));
   a.resize(k);
   for (int j = int(k - 1); j >= 0; j--)
   {
    a[j] *= b[0];
    for (int i = std::max(j - (int)n, 0) + 1; i < std::min(j + 1, (int)m); i++)
    {
     a[j] += a[j - i] * b[i];
    }
   }
  }
 }
 size_t com_size(size_t as, size_t bs)
 {
  if (!as || !bs)
  {
   return 0;
  }
  return std::max(flen, std::bit_ceil(as + bs - 1) / 2);
 }
 void mul_truncate(auto &a, auto const &b, size_t k)
 {
  using base = std::decay_t<decltype(a[0])>;
  if (std::min({k, std::size(a), std::size(b)}) < magic)
  {
   mul_slow(a, b, k);
   return;
  }
  auto n = std::max(flen, std::bit_ceil(
                              std::min(k, std::size(a)) + std::min(k, std::size(b)) - 1) /
                              2);
  auto A = dft<base>(a | std::views::take(k), n);
  auto B = dft<base>(b | std::views::take(k), n);
  a.resize((k + flen - 1) / flen * flen);
  A.mul_inplace(B, a, k);
  a.resize(k);
 }

 // store mod x^n-k in first half, x^n+k in second half
 void mod_split(auto &&x, size_t n, auto k)
 {
  using base = std::decay_t<decltype(k)>;
  dft<base>::init();
  assert(std::size(x) == 2 * n);
  u64x4 cur = u64x4{} + (k * bpow(base(2), 32)).getr();
  for (size_t i = 0; i < n; i += flen)
  {
   u64x4 xl = {
       x[i].getr(),
       x[i + 1].getr(),
       x[i + 2].getr(),
       x[i + 3].getr()};
   u64x4 xr = {
       x[n + i].getr(),
       x[n + i + 1].getr(),
       x[n + i + 2].getr(),
       x[n + i + 3].getr()};
   xr = montgomery_mul(xr, cur, dft<base>::mod, dft<base>::imod);
   xr = xr >= base::mod() ? xr - base::mod() : xr;
   auto t = xr;
   xr = xl - t;
   xl += t;
   xl = xl >= base::mod() ? xl - base::mod() : xl;
   xr = xr >= base::mod() ? xr + base::mod() : xr;
   for (size_t k = 0; k < flen; k++)
   {
    x[i + k].setr(typename base::UInt(xl[k]));
    x[n + i + k].setr(typename base::UInt(xr[k]));
   }
  }
  cp_algo::checkpoint("mod split");
 }
 void cyclic_mul(auto &a, auto &&b, size_t k)
 {
  assert(std::popcount(k) == 1);
  assert(std::size(a) == std::size(b) && std::size(a) == k);
  using base = std::decay_t<decltype(a[0])>;
  dft<base>::init();
  if (k <= (1 << 16))
  {
   auto ap = std::ranges::to<std::vector<base, big_alloc<base>>>(a);
   mul_truncate(ap, b, 2 * k);
   mod_split(ap, k, bpow(dft<base>::factor, k));
   std::ranges::copy(ap | std::views::take(k), begin(a));
   return;
  }
  k /= 2;
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
  for (size_t i = 0; i < k; i++)
  {
   auto t = (a[i] + a[i + k]) * i2;
   a[i + k] = (a[i] - a[i + k]) * factor;
   a[i] = t;
  }
  cp_algo::checkpoint("mod join");
 }
 auto make_copy(auto &&x)
 {
  return x;
 }
 void cyclic_mul(auto &a, auto const &b, size_t k)
 {
  return cyclic_mul(a, make_copy(b), k);
 }
 void mul(auto &a, auto &&b)
 {
  size_t N = size(a) + size(b);
  if (N > (1 << 20))
  {
   N--;
   size_t NN = std::bit_ceil(N);
   a.resize(NN);
   b.resize(NN);
   cyclic_mul(a, b, NN);
   a.resize(N);
  }
  else
  {
   mul_truncate(a, b, N - 1);
  }
 }
 void mul(auto &a, auto const &b)
 {
  size_t N = size(a) + size(b);
  if (N > (1 << 20))
  {
   mul(a, make_copy(b));
  }
  else
  {
   mul_truncate(a, b, N - 1);
  }
 }
}


// operations related to gcd and Euclidean algo
namespace cp_algo::math::poly::impl
{
 template <typename poly>
 using gcd_result = std::pair<
     std::list<std::decay_t<poly>>,
     linfrac<std::decay_t<poly>>>;

 template <typename poly>
 gcd_result<poly> half_gcd(poly &&A, poly &&B)
 {
  assert(A.deg() >= B.deg());
  size_t m = size(A.a) / 2;
  if (B.deg() < (int)m)
  {
   return {};
  }
  auto [ai, R] = A.divmod(B);
  std::tie(A, B) = {B, R};
  std::list a = {ai};
  auto T = -linfrac(ai).adj();

  auto advance = [&](size_t k)
  {
   auto [ak, Tk] = half_gcd(A.div_xk(k), B.div_xk(k));
   a.splice(end(a), ak);
   T.prepend(Tk);
   return Tk;
  };
  advance(m).apply(A, B);
  if constexpr (std::is_reference_v<poly>)
  {
   advance(2 * m - A.deg()).apply(A, B);
  }
  else
  {
   advance(2 * m - A.deg());
  }
  return {std::move(a), std::move(T)};
 }
 template <typename poly>
 gcd_result<poly> full_gcd(poly &&A, poly &&B)
 {
  using poly_t = std::decay_t<poly>;
  std::list<poly_t> ak;
  std::vector<linfrac<poly_t>> trs;
  while (!B.is_zero())
  {
   auto [a0, R] = A.divmod(B);
   ak.push_back(a0);
   trs.push_back(-linfrac(a0).adj());
   std::tie(A, B) = {B, R};

   auto [a, Tr] = half_gcd(A, B);
   ak.splice(end(ak), a);
   trs.push_back(Tr);
  }
  return {ak, std::accumulate(rbegin(trs), rend(trs), linfrac<poly_t>{}, std::multiplies{})};
 }

 // computes product of linfrac on [L, R)
 auto convergent(auto L, auto R)
 {
  using poly = decltype(L)::value_type;
  if (R == next(L))
  {
   return linfrac(*L);
  }
  else
  {
   int s = std::transform_reduce(L, R, 0, std::plus{}, std::mem_fn(&poly::deg));
   auto M = L;
   for (int c = M->deg(); 2 * c <= s; M++)
   {
    c += next(M)->deg();
   }
   return convergent(L, M) * convergent(M, R);
  }
 }
 template <typename poly>
 poly min_rec(poly const &p, size_t d)
 {
  auto R2 = p.mod_xk(d).reversed(d), R1 = poly::xk(d);
  if (R2.is_zero())
  {
   return poly(1);
  }
  auto [a, Tr] = full_gcd(R1, R2);
  a.emplace_back();
  auto pref = begin(a);
  for (int delta = (int)d - a.front().deg(); delta >= 0; pref++)
  {
   delta -= pref->deg() + next(pref)->deg();
  }
  return convergent(begin(a), pref).a;
 }

 template <typename poly>
 std::optional<poly> inv_mod(poly p, poly q)
 {
  assert(!q.is_zero());
  auto [a, Tr] = full_gcd(q, p);
  if (q.deg() != 0)
  {
   return std::nullopt;
  }
  return Tr.b / q[0];
 }
}




// operations related to polynomial division
namespace cp_algo::math::poly::impl
{
 auto divmod_slow(auto const &p, auto const &q)
 {
  auto R = p;
  auto D = decltype(p){};
  auto q_lead_inv = q.lead().inv();
  while (R.deg() >= q.deg())
  {
   D.a.push_back(R.lead() * q_lead_inv);
   if (D.lead() != 0)
   {
    for (size_t i = 1; i <= q.a.size(); i++)
    {
     R.a[R.a.size() - i] -= D.lead() * q.a[q.a.size() - i];
    }
   }
   R.a.pop_back();
  }
  std::ranges::reverse(D.a);
  R.normalize();
  return std::array{D, R};
 }
 template <typename poly>
 auto divmod_hint(poly const &p, poly const &q, poly const &qri)
 {
  assert(!q.is_zero());
  int d = p.deg() - q.deg();
  if (std::min(d, q.deg()) < magic)
  {
   return divmod_slow(p, q);
  }
  poly D;
  if (d >= 0)
  {
   D = (p.reversed().mod_xk(d + 1) * qri.mod_xk(d + 1)).mod_xk(d + 1).reversed(d + 1);
  }
  return std::array{D, p - D * q};
 }
 auto divmod(auto const &p, auto const &q)
 {
  assert(!q.is_zero());
  int d = p.deg() - q.deg();
  if (std::min(d, q.deg()) < magic)
  {
   return divmod_slow(p, q);
  }
  return divmod_hint(p, q, q.reversed().inv(d + 1));
 }

 template <typename poly>
 poly powmod_hint(poly const &p, int64_t k, poly const &md, poly const &mdri)
 {
  return bpow(p, k, poly(1), [&](auto const &p, auto const &q)
              { return divmod_hint(p * q, md, mdri)[1]; });
 }
 template <typename poly>
 auto powmod(poly const &p, int64_t k, poly const &md)
 {
  int d = md.deg();
  if (p == poly::xk(1) && false)
  { // does it actually speed anything up?..
   if (k < md.deg())
   {
    return poly::xk(k);
   }
   else
   {
    auto mdr = md.reversed();
    return (mdr.inv(k - md.deg() + 1, md.deg()) * mdr).reversed(md.deg());
   }
  }
  if (md == poly::xk(d))
  {
   return p.pow(k, d);
  }
  if (md == poly::xk(d) - poly(1))
  {
   return p.powmod_circular(k, d);
  }
  return powmod_hint(p, k, md, md.reversed().inv(md.deg() + 1));
 }
 template <typename poly>
 poly &inv_inplace(poly &q, int64_t k, size_t n)
 {
  using poly_t = std::decay_t<poly>;
  using base = poly_t::base;
  if (k <= std::max<int64_t>(n, size(q.a)))
  {
   return q.inv_inplace(k + n).div_xk_inplace(k);
  }
  if (k % 2)
  {
   return inv_inplace(q, k - 1, n + 1).div_xk_inplace(1);
  }
  auto [q0, q1] = q.bisect();
  auto qq = q0 * q0 - (q1 * q1).mul_xk_inplace(1);
  inv_inplace(qq, k / 2 - q.deg() / 2, (n + 1) / 2 + q.deg() / 2);
  size_t N = fft::com_size(size(q0.a), size(qq.a));
  auto q0f = fft::dft<base>(q0.a, N);
  auto q1f = fft::dft<base>(q1.a, N);
  auto qqf = fft::dft<base>(qq.a, N);
  size_t M = q0.deg() + (n + 1) / 2;
  typename poly::Vector A, B;
  A.resize((M + fft::flen - 1) / fft::flen * fft::flen);
  B.resize((M + fft::flen - 1) / fft::flen * fft::flen);
  q0f.mul(qqf, A, M);
  q1f.mul_inplace(qqf, B, M);
  q.a.resize(n + 1);
  for (size_t i = 0; i < n; i += 2)
  {
   q.a[i] = A[q0.deg() + i / 2];
   q.a[i + 1] = -B[q0.deg() + i / 2];
  }
  q.a.pop_back();
  q.normalize();
  return q;
 }
 template <typename poly>
 poly &inv_inplace(poly &p, size_t n)
 {
  using poly_t = std::decay_t<poly>;
  using base = poly_t::base;
  if (n == 1)
  {
   return p = base(1) / p[0];
  }
  // Q(-x) = P0(x^2) + xP1(x^2)
  auto [q0, q1] = p.bisect(n);

  size_t N = fft::com_size(size(q0.a), (n + 1) / 2);

  auto q0f = fft::dft<base>(q0.a, N);
  auto q1f = fft::dft<base>(q1.a, N);

  // Q(x)*Q(-x) = Q0(x^2)^2 - x^2 Q1(x^2)^2
  auto qq = poly_t(q0f * q0f) - poly_t(q1f * q1f).mul_xk_inplace(1);

  inv_inplace(qq, (n + 1) / 2);
  auto qqf = fft::dft<base>(qq.a, N);

  typename poly::Vector A, B;
  A.resize(((n + 1) / 2 + fft::flen - 1) / fft::flen * fft::flen);
  B.resize(((n + 1) / 2 + fft::flen - 1) / fft::flen * fft::flen);
  q0f.mul(qqf, A, (n + 1) / 2);
  q1f.mul_inplace(qqf, B, (n + 1) / 2);
  p.a.resize(n + 1);
  for (size_t i = 0; i < n; i += 2)
  {
   p.a[i] = A[i / 2];
   p.a[i + 1] = -B[i / 2];
  }
  p.a.pop_back();
  p.normalize();
  return p;
 }
}




namespace cp_algo::math
{
 // fact/rfact/small_inv are caching
 // Beware of usage with dynamic mod
 template <typename T>
 T fact(int n)
 {
  static std::vector<T> F(maxn);
  static bool init = false;
  if (!init)
  {
   F[0] = T(1);
   for (int i = 1; i < maxn; i++)
   {
    F[i] = F[i - 1] * T(i);
   }
   init = true;
  }
  return F[n];
 }
 // Only works for modint types
 template <typename T>
 T rfact(int n)
 {
  static std::vector<T> F(maxn);
  static bool init = false;
  if (!init)
  {
   int t = (int)std::min<int64_t>(T::mod(), maxn) - 1;
   F[t] = T(1) / fact<T>(t);
   for (int i = t - 1; i >= 0; i--)
   {
    F[i] = F[i + 1] * T(i + 1);
   }
   init = true;
  }
  return F[n];
 }
 template <typename T>
 T small_inv(int n)
 {
  static std::vector<T> F(maxn);
  static bool init = false;
  if (!init)
  {
   for (int i = 1; i < maxn; i++)
   {
    F[i] = rfact<T>(i) * fact<T>(i - 1);
   }
   init = true;
  }
  return F[n];
 }
 template <typename T>
 T binom_large(T n, int r)
 {
  assert(r < maxn);
  T ans = 1;
  for (int i = 0; i < r; i++)
  {
   ans = ans * T(n - i) * small_inv<T>(i + 1);
  }
  return ans;
 }
 template <typename T>
 T binom(int n, int r)
 {
  if (r < 0 || r > n)
  {
   return T(0);
  }
  else if (n >= maxn)
  {
   return binom_large(T(n), r);
  }
  else
  {
   return fact<T>(n) * rfact<T>(r) * rfact<T>(n - r);
  }
 }
}




namespace cp_algo::math
{
 // https://en.wikipedia.org/wiki/Berlekamp-Rabin_algorithm
 template <modint_type base>
 std::optional<base> sqrt(base b)
 {
  if (b == base(0))
  {
   return base(0);
  }
  else if (bpow(b, (b.mod() - 1) / 2) != base(1))
  {
   return std::nullopt;
  }
  else
  {
   while (true)
   {
    base z = random::rng();
    if (z * z == b)
    {
     return z;
    }
    lin<base> x(1, z, b); // x + z (mod x^2 - b)
    x = bpow(x, (b.mod() - 1) / 2, lin<base>(0, 1, b));
    if (x.a != base(0))
    {
     return x.a.inv();
    }
   }
  }
 }
}


namespace cp_algo::math
{
 template <typename T, class Alloc = big_alloc<T>>
 struct poly_t
 {
  using Vector = std::vector<T, Alloc>;
  using base = T;
  Vector a;

  poly_t &normalize()
  {
   while (deg() >= 0 && lead() == base(0))
   {
    a.pop_back();
   }
   return *this;
  }

  poly_t() {}
  poly_t(T a0) : a{a0} { normalize(); }
  poly_t(Vector const &t) : a(t) { normalize(); }
  poly_t(Vector &&t) : a(std::move(t)) { normalize(); }

  poly_t &negate_inplace()
  {
   std::ranges::transform(a, begin(a), std::negate{});
   return *this;
  }
  poly_t operator-() const
  {
   return poly_t(*this).negate_inplace();
  }
  poly_t &operator+=(poly_t const &t)
  {
   a.resize(std::max(size(a), size(t.a)));
   std::ranges::transform(a, t.a, begin(a), std::plus{});
   return normalize();
  }
  poly_t &operator-=(poly_t const &t)
  {
   a.resize(std::max(size(a), size(t.a)));
   std::ranges::transform(a, t.a, begin(a), std::minus{});
   return normalize();
  }
  poly_t operator+(poly_t const &t) const { return poly_t(*this) += t; }
  poly_t operator-(poly_t const &t) const { return poly_t(*this) -= t; }

  poly_t &mod_xk_inplace(size_t k)
  {
   a.resize(std::min(size(a), k));
   return normalize();
  }
  poly_t &mul_xk_inplace(size_t k)
  {
   a.insert(begin(a), k, T(0));
   return normalize();
  }
  poly_t &div_xk_inplace(int64_t k)
  {
   if (k < 0)
   {
    return mul_xk_inplace(-k);
   }
   a.erase(begin(a), begin(a) + std::min<size_t>(k, size(a)));
   return normalize();
  }
  poly_t &substr_inplace(size_t l, size_t k)
  {
   return mod_xk_inplace(l + k).div_xk_inplace(l);
  }
  poly_t mod_xk(size_t k) const { return poly_t(*this).mod_xk_inplace(k); }
  poly_t mul_xk(size_t k) const { return poly_t(*this).mul_xk_inplace(k); }
  poly_t div_xk(int64_t k) const { return poly_t(*this).div_xk_inplace(k); }
  poly_t substr(size_t l, size_t k) const { return poly_t(*this).substr_inplace(l, k); }

  poly_t &operator*=(const poly_t &t)
  {
   fft::mul(a, t.a);
   normalize();
   return *this;
  }
  poly_t operator*(const poly_t &t) const { return poly_t(*this) *= t; }

  poly_t &operator/=(const poly_t &t) { return *this = divmod(t)[0]; }
  poly_t &operator%=(const poly_t &t) { return *this = divmod(t)[1]; }
  poly_t operator/(poly_t const &t) const { return poly_t(*this) /= t; }
  poly_t operator%(poly_t const &t) const { return poly_t(*this) %= t; }

  poly_t &operator*=(T const &x)
  {
   for (auto &it : a)
   {
    it *= x;
   }
   return normalize();
  }
  poly_t &operator/=(T const &x) { return *this *= x.inv(); }
  poly_t operator*(T const &x) const { return poly_t(*this) *= x; }
  poly_t operator/(T const &x) const { return poly_t(*this) /= x; }

  poly_t &reverse(size_t n)
  {
   a.resize(n);
   std::ranges::reverse(a);
   return normalize();
  }
  poly_t &reverse() { return reverse(size(a)); }
  poly_t reversed(size_t n) const { return poly_t(*this).reverse(n); }
  poly_t reversed() const { return poly_t(*this).reverse(); }

  std::array<poly_t, 2> divmod(poly_t const &b) const
  {
   return poly::impl::divmod(*this, b);
  }

  // reduces A/B to A'/B' such that
  // deg B' < deg A / 2
  static std::pair<std::list<poly_t>, linfrac<poly_t>> half_gcd(auto &&A, auto &&B)
  {
   return poly::impl::half_gcd(A, B);
  }
  // reduces A / B to gcd(A, B) / 0
  static std::pair<std::list<poly_t>, linfrac<poly_t>> full_gcd(auto &&A, auto &&B)
  {
   return poly::impl::full_gcd(A, B);
  }
  static poly_t gcd(poly_t &&A, poly_t &&B)
  {
   full_gcd(A, B);
   return A;
  }

  // Returns a (non-monic) characteristic polynomial
  // of the minimum linear recurrence for the sequence
  poly_t min_rec(size_t d) const
  {
   return poly::impl::min_rec(*this, d);
  }

  // calculate inv to *this modulo t
  std::optional<poly_t> inv_mod(poly_t const &t) const
  {
   return poly::impl::inv_mod(*this, t);
  };

  poly_t negx() const
  { // A(x) -> A(-x)
   auto res = *this;
   for (int i = 1; i <= deg(); i += 2)
   {
    res.a[i] = -res[i];
   }
   return res;
  }

  void print(int n) const
  {
   for (int i = 0; i < n; i++)
   {
    std::cout << (*this)[i] << ' ';
   }
   std::cout << "\n";
  }

  void print() const
  {
   print(deg() + 1);
  }

  T eval(T x) const
  { // evaluates in single point x
   T res(0);
   for (int i = deg(); i >= 0; i--)
   {
    res *= x;
    res += a[i];
   }
   return res;
  }

  T lead() const
  { // leading coefficient
   assert(!is_zero());
   return a.back();
  }

  int deg() const
  { // degree, -1 for P(x) = 0
   return (int)a.size() - 1;
  }

  bool is_zero() const
  {
   return a.empty();
  }

  T operator[](int idx) const
  {
   return idx < 0 || idx > deg() ? T(0) : a[idx];
  }

  T &coef(size_t idx)
  { // mutable reference at coefficient
   return a[idx];
  }

  bool operator==(const poly_t &t) const { return a == t.a; }
  bool operator!=(const poly_t &t) const { return a != t.a; }

  poly_t &deriv_inplace(int k = 1)
  {
   if (deg() + 1 < k)
   {
    return *this = poly_t{};
   }
   for (int i = k; i <= deg(); i++)
   {
    a[i - k] = fact<T>(i) * rfact<T>(i - k) * a[i];
   }
   a.resize(deg() + 1 - k);
   return *this;
  }
  poly_t deriv(int k = 1) const
  { // calculate derivative
   return poly_t(*this).deriv_inplace(k);
  }

  poly_t &integr_inplace()
  {
   a.push_back(0);
   for (int i = deg() - 1; i >= 0; i--)
   {
    a[i + 1] = a[i] * small_inv<T>(i + 1);
   }
   a[0] = 0;
   return *this;
  }
  poly_t integr() const
  { // calculate integral with C = 0
   Vector res(deg() + 2);
   for (int i = 0; i <= deg(); i++)
   {
    res[i + 1] = a[i] * small_inv<T>(i + 1);
   }
   return res;
  }

  size_t trailing_xk() const
  { // Let p(x) = x^k * t(x), return k
   if (is_zero())
   {
    return -1;
   }
   int res = 0;
   while (a[res] == T(0))
   {
    res++;
   }
   return res;
  }

  // calculate log p(x) mod x^n
  poly_t &log_inplace(size_t n)
  {
   assert(a[0] == T(1));
   mod_xk_inplace(n);
   return (inv_inplace(n) *= mod_xk_inplace(n).deriv()).mod_xk_inplace(n - 1).integr_inplace();
  }
  poly_t log(size_t n) const
  {
   return poly_t(*this).log_inplace(n);
  }

  poly_t &mul_truncate(poly_t const &t, size_t k)
  {
   fft::mul_truncate(a, t.a, k);
   return normalize();
  }

  poly_t &exp_inplace(size_t n)
  {
   if (is_zero())
   {
    return *this = T(1);
   }
   assert(a[0] == T(0));
   a[0] = 1;
   size_t a = 1;
   while (a < n)
   {
    poly_t C = log(2 * a).div_xk_inplace(a) - substr(a, 2 * a);
    *this -= C.mul_truncate(*this, a).mul_xk_inplace(a);
    a *= 2;
   }
   return mod_xk_inplace(n);
  }

  poly_t exp(size_t n) const
  { // calculate exp p(x) mod x^n
   return poly_t(*this).exp_inplace(n);
  }

  poly_t pow_bin(int64_t k, size_t n) const
  { // O(n log n log k)
   if (k == 0)
   {
    return poly_t(1).mod_xk(n);
   }
   else
   {
    auto t = pow(k / 2, n);
    t = (t * t).mod_xk(n);
    return (k % 2 ? *this * t : t).mod_xk(n);
   }
  }

  poly_t circular_closure(size_t m) const
  {
   if (deg() == -1)
   {
    return *this;
   }
   auto t = *this;
   for (size_t i = t.deg(); i >= m; i--)
   {
    t.a[i - m] += t.a[i];
   }
   t.a.resize(std::min(t.a.size(), m));
   return t;
  }

  static poly_t mul_circular(poly_t const &a, poly_t const &b, size_t m)
  {
   return (a.circular_closure(m) * b.circular_closure(m)).circular_closure(m);
  }

  poly_t powmod_circular(int64_t k, size_t m) const
  {
   if (k == 0)
   {
    return poly_t(1);
   }
   else
   {
    auto t = powmod_circular(k / 2, m);
    t = mul_circular(t, t, m);
    if (k % 2)
    {
     t = mul_circular(t, *this, m);
    }
    return t;
   }
  }

  poly_t powmod(int64_t k, poly_t const &md) const
  {
   return poly::impl::powmod(*this, k, md);
  }

  // O(d * n) with the derivative trick from
  // https://codeforces.com/blog/entry/73947?#comment-581173
  poly_t pow_dn(int64_t k, size_t n) const
  {
   if (n == 0)
   {
    return poly_t(T(0));
   }
   assert((*this)[0] != T(0));
   Vector Q(n);
   Q[0] = bpow(a[0], k);
   auto a0inv = a[0].inv();
   for (int i = 1; i < (int)n; i++)
   {
    for (int j = 1; j <= std::min(deg(), i); j++)
    {
     Q[i] += a[j] * Q[i - j] * (T(k) * T(j) - T(i - j));
    }
    Q[i] *= small_inv<T>(i) * a0inv;
   }
   return Q;
  }

  // calculate p^k(n) mod x^n in O(n log n)
  // might be quite slow due to high constant
  poly_t pow(int64_t k, size_t n) const
  {
   if (is_zero())
   {
    return k ? *this : poly_t(1);
   }
   size_t i = trailing_xk();
   if (i > 0)
   {
    return k >= int64_t(n + i - 1) / (int64_t)i ? poly_t(T(0)) : div_xk(i).pow(k, n - i * k).mul_xk(i * k);
   }
   if (std::min(deg(), (int)n) <= magic)
   {
    return pow_dn(k, n);
   }
   if (k <= magic)
   {
    return pow_bin(k, n);
   }
   T j = a[i];
   poly_t t = *this / j;
   return bpow(j, k) * (t.log(n) * T(k)).exp(n).mod_xk(n);
  }

  // returns std::nullopt if undefined
  std::optional<poly_t> sqrt(size_t n) const
  {
   if (is_zero())
   {
    return *this;
   }
   size_t i = trailing_xk();
   if (i % 2)
   {
    return std::nullopt;
   }
   else if (i > 0)
   {
    auto ans = div_xk(i).sqrt(n - i / 2);
    return ans ? ans->mul_xk(i / 2) : ans;
   }
   auto st = math::sqrt((*this)[0]);
   if (st)
   {
    poly_t ans = *st;
    size_t a = 1;
    while (a < n)
    {
     a *= 2;
     ans -= (ans - mod_xk(a) * ans.inv(a)).mod_xk(a) / 2;
    }
    return ans.mod_xk(n);
   }
   return std::nullopt;
  }

  poly_t mulx(T a) const
  { // component-wise multiplication with a^k
   T cur = 1;
   poly_t res(*this);
   for (int i = 0; i <= deg(); i++)
   {
    res.coef(i) *= cur;
    cur *= a;
   }
   return res;
  }

  poly_t mulx_sq(T a) const
  { // component-wise multiplication with a^{k choose 2}
   T cur = 1, total = 1;
   poly_t res(*this);
   for (int i = 0; i <= deg(); i++)
   {
    res.coef(i) *= total;
    cur *= a;
    total *= cur;
   }
   return res;
  }

  // be mindful of maxn, as the function
  // requires multiplying polynomials of size deg() and n+deg()!
  poly_t chirpz(T z, int n) const
  { // P(1), P(z), P(z^2), ..., P(z^(n-1))
   if (is_zero())
   {
    return Vector(n, 0);
   }
   if (z == T(0))
   {
    Vector ans(n, (*this)[0]);
    if (n > 0)
    {
     ans[0] = accumulate(begin(a), end(a), T(0));
    }
    return ans;
   }
   auto A = mulx_sq(z.inv());
   auto B = ones(n + deg()).mulx_sq(z);
   return semicorr(B, A).mod_xk(n).mulx_sq(z.inv());
  }

  // res[i] = prod_{1 <= j <= i} 1/(1 - z^j)
  static auto _1mzk_prod_inv(T z, int n)
  {
   Vector res(n, 1), zk(n);
   zk[0] = 1;
   for (int i = 1; i < n; i++)
   {
    zk[i] = zk[i - 1] * z;
    res[i] = res[i - 1] * (T(1) - zk[i]);
   }
   res.back() = res.back().inv();
   for (int i = n - 2; i >= 0; i--)
   {
    res[i] = (T(1) - zk[i + 1]) * res[i + 1];
   }
   return res;
  }

  // prod_{0 <= j < n} (1 - z^j x)
  static auto _1mzkx_prod(T z, int n)
  {
   if (n == 1)
   {
    return poly_t(Vector{1, -1});
   }
   else
   {
    auto t = _1mzkx_prod(z, n / 2);
    t *= t.mulx(bpow(z, n / 2));
    if (n % 2)
    {
     t *= poly_t(Vector{1, -bpow(z, n - 1)});
    }
    return t;
   }
  }

  poly_t chirpz_inverse(T z, int n) const
  { // P(1), P(z), P(z^2), ..., P(z^(n-1))
   if (is_zero())
   {
    return {};
   }
   if (z == T(0))
   {
    if (n == 1)
    {
     return *this;
    }
    else
    {
     return Vector{(*this)[1], (*this)[0] - (*this)[1]};
    }
   }
   Vector y(n);
   for (int i = 0; i < n; i++)
   {
    y[i] = (*this)[i];
   }
   auto prods_pos = _1mzk_prod_inv(z, n);
   auto prods_neg = _1mzk_prod_inv(z.inv(), n);

   T zn = bpow(z, n - 1).inv();
   T znk = 1;
   for (int i = 0; i < n; i++)
   {
    y[i] *= znk * prods_neg[i] * prods_pos[(n - 1) - i];
    znk *= zn;
   }

   poly_t p_over_q = poly_t(y).chirpz(z, n);
   poly_t q = _1mzkx_prod(z, n);

   return (p_over_q * q).mod_xk_inplace(n).reverse(n);
  }

  static poly_t build(std::vector<poly_t> &res, int v, auto L, auto R)
  { // builds evaluation tree for (x-a1)(x-a2)...(x-an)
   if (R - L == 1)
   {
    return res[v] = Vector{-*L, 1};
   }
   else
   {
    auto M = L + (R - L) / 2;
    return res[v] = build(res, 2 * v, L, M) * build(res, 2 * v + 1, M, R);
   }
  }

  poly_t to_newton(std::vector<poly_t> &tree, int v, auto l, auto r)
  {
   if (r - l == 1)
   {
    return *this;
   }
   else
   {
    auto m = l + (r - l) / 2;
    auto A = (*this % tree[2 * v]).to_newton(tree, 2 * v, l, m);
    auto B = (*this / tree[2 * v]).to_newton(tree, 2 * v + 1, m, r);
    return A + B.mul_xk(m - l);
   }
  }

  poly_t to_newton(Vector p)
  {
   if (is_zero())
   {
    return *this;
   }
   size_t n = p.size();
   std::vector<poly_t> tree(4 * n);
   build(tree, 1, begin(p), end(p));
   return to_newton(tree, 1, begin(p), end(p));
  }

  Vector eval(std::vector<poly_t> &tree, int v, auto l, auto r)
  { // auxiliary evaluation function
   if (r - l == 1)
   {
    return {eval(*l)};
   }
   else
   {
    auto m = l + (r - l) / 2;
    auto A = (*this % tree[2 * v]).eval(tree, 2 * v, l, m);
    auto B = (*this % tree[2 * v + 1]).eval(tree, 2 * v + 1, m, r);
    A.insert(end(A), begin(B), end(B));
    return A;
   }
  }

  Vector eval(Vector x)
  { // evaluate polynomial in (x1, ..., xn)
   size_t n = x.size();
   if (is_zero())
   {
    return Vector(n, T(0));
   }
   std::vector<poly_t> tree(4 * n);
   build(tree, 1, begin(x), end(x));
   return eval(tree, 1, begin(x), end(x));
  }

  poly_t inter(std::vector<poly_t> &tree, int v, auto ly, auto ry)
  { // auxiliary interpolation function
   if (ry - ly == 1)
   {
    return {*ly / a[0]};
   }
   else
   {
    auto my = ly + (ry - ly) / 2;
    auto A = (*this % tree[2 * v]).inter(tree, 2 * v, ly, my);
    auto B = (*this % tree[2 * v + 1]).inter(tree, 2 * v + 1, my, ry);
    return A * tree[2 * v + 1] + B * tree[2 * v];
   }
  }

  static auto inter(Vector x, Vector y)
  { // interpolates minimum polynomial from (xi, yi) pairs
   size_t n = x.size();
   std::vector<poly_t> tree(4 * n);
   return build(tree, 1, begin(x), end(x)).deriv().inter(tree, 1, begin(y), end(y));
  }

  static auto resultant(poly_t a, poly_t b)
  { // computes resultant of a and b
   if (b.is_zero())
   {
    return 0;
   }
   else if (b.deg() == 0)
   {
    return bpow(b.lead(), a.deg());
   }
   else
   {
    int pw = a.deg();
    a %= b;
    pw -= a.deg();
    auto mul = bpow(b.lead(), pw) * T((b.deg() & a.deg() & 1) ? -1 : 1);
    auto ans = resultant(b, a);
    return ans * mul;
   }
  }

  static poly_t xk(size_t n)
  { // P(x) = x^n
   return poly_t(T(1)).mul_xk(n);
  }

  static poly_t ones(size_t n)
  { // P(x) = 1 + x + ... + x^{n-1}
   return Vector(n, 1);
  }

  static poly_t expx(size_t n)
  { // P(x) = e^x (mod x^n)
   return ones(n).borel();
  }

  static poly_t log1px(size_t n)
  { // P(x) = log(1+x) (mod x^n)
   Vector coeffs(n, 0);
   for (size_t i = 1; i < n; i++)
   {
    coeffs[i] = (i & 1 ? T(i).inv() : -T(i).inv());
   }
   return coeffs;
  }

  static poly_t log1mx(size_t n)
  { // P(x) = log(1-x) (mod x^n)
   return -ones(n).integr();
  }

  // [x^k] (a corr b) = sum_{i} a{(k-m)+i}*bi
  static poly_t corr(poly_t const &a, poly_t const &b)
  { // cross-correlation
   return a * b.reversed();
  }

  // [x^k] (a semicorr b) = sum_i a{i+k} * b{i}
  static poly_t semicorr(poly_t const &a, poly_t const &b)
  {
   return corr(a, b).div_xk(b.deg());
  }

  poly_t invborel() const
  { // ak *= k!
   auto res = *this;
   for (int i = 0; i <= deg(); i++)
   {
    res.coef(i) *= fact<T>(i);
   }
   return res;
  }

  poly_t borel() const
  { // ak /= k!
   auto res = *this;
   for (int i = 0; i <= deg(); i++)
   {
    res.coef(i) *= rfact<T>(i);
   }
   return res;
  }

  poly_t shift(T a) const
  { // P(x + a)
   return semicorr(invborel(), expx(deg() + 1).mulx(a)).borel();
  }

  poly_t x2()
  { // P(x) -> P(x^2)
   Vector res(2 * a.size());
   for (size_t i = 0; i < a.size(); i++)
   {
    res[2 * i] = a[i];
   }
   return res;
  }

  // Return {P0, P1}, where P(x) = P0(x) + xP1(x)
  std::array<poly_t, 2> bisect(size_t n) const
  {
   n = std::min(n, size(a));
   Vector res[2];
   for (size_t i = 0; i < n; i++)
   {
    res[i % 2].push_back(a[i]);
   }
   return {res[0], res[1]};
  }
  std::array<poly_t, 2> bisect() const
  {
   return bisect(size(a));
  }

  // Find [x^k] P / Q
  static T kth_rec_inplace(poly_t &P, poly_t &Q, int64_t k)
  {
   while (k > Q.deg())
   {
    size_t n = Q.a.size();
    auto [Q0, Q1] = Q.bisect();
    auto [P0, P1] = P.bisect();

    size_t N = fft::com_size((n + 1) / 2, (n + 1) / 2);

    auto Q0f = fft::dft<T>(Q0.a, N);
    auto Q1f = fft::dft<T>(Q1.a, N);
    auto P0f = fft::dft<T>(P0.a, N);
    auto P1f = fft::dft<T>(P1.a, N);

    Q = poly_t(Q0f * Q0f) -= poly_t(Q1f * Q1f).mul_xk_inplace(1);
    if (k % 2)
    {
     P = poly_t(Q0f *= P1f) -= poly_t(Q1f *= P0f);
    }
    else
    {
     P = poly_t(Q0f *= P0f) -= poly_t(Q1f *= P1f).mul_xk_inplace(1);
    }
    k /= 2;
   }
   return (P *= Q.inv_inplace(Q.deg() + 1))[(int)k];
  }
  static T kth_rec(poly_t const &P, poly_t const &Q, int64_t k)
  {
   return kth_rec_inplace(poly_t(P), poly_t(Q), k);
  }

  // inverse series mod x^n
  poly_t &inv_inplace(size_t n)
  {
   return poly::impl::inv_inplace(*this, n);
  }
  poly_t inv(size_t n) const
  {
   return poly_t(*this).inv_inplace(n);
  }
  // [x^k]..[x^{k+n-1}] of inv()
  // supports negative k if k+n >= 0
  poly_t &inv_inplace(int64_t k, size_t n)
  {
   return poly::impl::inv_inplace(*this, k, n);
  }
  poly_t inv(int64_t k, size_t n) const
  {
   return poly_t(*this).inv_inplace(k, n);
  }

  // compute A(B(x)) mod x^n in O(n^2)
  static poly_t compose(poly_t A, poly_t B, int n)
  {
   int q = std::sqrt(n);
   std::vector<poly_t> Bk(q);
   auto Bq = B.pow(q, n);
   Bk[0] = poly_t(T(1));
   for (int i = 1; i < q; i++)
   {
    Bk[i] = (Bk[i - 1] * B).mod_xk(n);
   }
   poly_t Bqk(1);
   poly_t ans;
   for (int i = 0; i <= n / q; i++)
   {
    poly_t cur;
    for (int j = 0; j < q; j++)
    {
     cur += Bk[j] * A[i * q + j];
    }
    ans += (Bqk * cur).mod_xk(n);
    Bqk = (Bqk * Bq).mod_xk(n);
   }
   return ans;
  }

  // compute A(B(x)) mod x^n in O(sqrt(pqn log^3 n))
  // preferrable when p = deg A and q = deg B
  // are much less than n
  static poly_t compose_large(poly_t A, poly_t B, int n)
  {
   if (B[0] != T(0))
   {
    return compose_large(A.shift(B[0]), B - B[0], n);
   }

   int q = std::sqrt(n);
   auto [B0, B1] = std::make_pair(B.mod_xk(q), B.div_xk(q));

   B0 = B0.div_xk(1);
   std::vector<poly_t> pw(A.deg() + 1);
   auto getpow = [&](int k)
   {
    return pw[k].is_zero() ? pw[k] = B0.pow(k, n - k) : pw[k];
   };

   std::function<poly_t(poly_t const &, int, int)> compose_dac = [&getpow, &compose_dac](poly_t const &f, int m, int N)
   {
    if (f.deg() <= 0)
    {
     return f;
    }
    int k = m / 2;
    auto [f0, f1] = std::make_pair(f.mod_xk(k), f.div_xk(k));
    auto [A, B] = std::make_pair(compose_dac(f0, k, N), compose_dac(f1, m - k, N - k));
    return (A + (B.mod_xk(N - k) * getpow(k).mod_xk(N - k)).mul_xk(k)).mod_xk(N);
   };

   int r = n / q;
   auto Ar = A.deriv(r);
   auto AB0 = compose_dac(Ar, Ar.deg() + 1, n);

   auto Bd = B0.mul_xk(1).deriv();

   poly_t ans = T(0);

   std::vector<poly_t> B1p(r + 1);
   B1p[0] = poly_t(T(1));
   for (int i = 1; i <= r; i++)
   {
    B1p[i] = (B1p[i - 1] * B1.mod_xk(n - i * q)).mod_xk(n - i * q);
   }
   while (r >= 0)
   {
    ans += (AB0.mod_xk(n - r * q) * rfact<T>(r) * B1p[r]).mul_xk(r * q).mod_xk(n);
    r--;
    if (r >= 0)
    {
     AB0 = ((AB0 * Bd).integr() + A[r] * fact<T>(r)).mod_xk(n);
    }
   }

   return ans;
  }
 };
 template <typename base>
 static auto operator*(const auto &a, const poly_t<base> &b)
 {
  return b * a;
 }
};






namespace cp_algo::linalg
{
 template <typename base, class Alloc = big_alloc<base>>
 struct vec : std::basic_string<base, std::char_traits<base>, Alloc>
 {
  using Base = std::basic_string<base, std::char_traits<base>, Alloc>;
  using Base::Base;

  vec(Base const &t) : Base(t) {}
  vec(Base &&t) : Base(std::move(t)) {}
  vec(size_t n) : Base(n, base()) {}
  vec(auto &&r) : Base(std::ranges::to<Base>(r)) {}

  static vec ei(size_t n, size_t i)
  {
   vec res(n);
   res[i] = 1;
   return res;
  }

  auto operator-() const
  {
   return *this | std::views::transform([](auto x)
                                        { return -x; });
  }
  auto operator*(base t) const
  {
   return *this | std::views::transform([t](auto x)
                                        { return x * t; });
  }
  auto operator*=(base t)
  {
   for (auto &it : *this)
   {
    it *= t;
   }
   return *this;
  }

  virtual void add_scaled(vec const &b, base scale, size_t i = 0)
  {
   if (scale != base(0))
   {
    for (; i < size(*this); i++)
    {
     (*this)[i] += scale * b[i];
    }
   }
  }
  virtual vec const &normalize()
  {
   return static_cast<vec &>(*this);
  }
  virtual base normalize(size_t i)
  {
   return (*this)[i];
  }
  void read()
  {
   for (auto &it : *this)
   {
    std::cin >> it;
   }
  }
  void print() const
  {
   for (auto &it : *this)
   {
    std::cout << it << " ";
   }
   std::cout << "\n";
  }
  static vec random(size_t n)
  {
   vec res(n);
   std::ranges::generate(res, random::rng);
   return res;
  }
  // Concatenate vectors
  vec operator|(vec const &t) const
  {
   return std::views::join(std::array{
       std::views::all(*this),
       std::views::all(t)});
  }

  // Generally, vec shouldn't be modified
  // after its pivot index is set
  std::pair<size_t, base> find_pivot()
  {
   if (pivot == size_t(-1))
   {
    pivot = 0;
    while (pivot < size(*this) && normalize(pivot) == base(0))
    {
     pivot++;
    }
    if (pivot < size(*this))
    {
     pivot_inv = base(1) / (*this)[pivot];
    }
   }
   return {pivot, pivot_inv};
  }
  void reduce_by(vec &t)
  {
   auto [pivot, pinv] = t.find_pivot();
   if (pivot < size(*this))
   {
    add_scaled(t, -normalize(pivot) * pinv, pivot);
   }
  }

 private:
  size_t pivot = -1;
  base pivot_inv;
 };

 template <math::modint_type base, class Alloc = big_alloc<base>>
 struct modint_vec : vec<base, Alloc>
 {
  using Base = vec<base, Alloc>;
  using Base::Base;

  modint_vec(Base const &t) : Base(t) {}
  modint_vec(Base &&t) : Base(std::move(t)) {}

  void add_scaled(Base const &b, base scale, size_t i = 0) override
  {
   static_assert(base::bits >= 64, "Only wide modint types for linalg");
   if (scale != base(0))
   {
    assert(Base::size() == b.size());
    size_t n = size(*this);
    u64x4 scaler = u64x4() + scale.getr();
    if (is_aligned(&(*this)[0]) && is_aligned(&b[0])) // verify we're not in SSO
     for (i -= i % 4; i < n; i += 4)
     {
      auto &ai = vector_cast<u64x4>((*this)[i]);
      auto bi = vector_cast<u64x4 const>(b[i]);
#ifdef __AVX2__
      ai += u64x4(_mm256_mul_epu32(__m256i(scaler), __m256i(bi)));
#else
      ai += scaler * bi;
#endif
     }
    for (; i < n; i++)
    {
     (*this)[i].add_unsafe(b[i].getr_direct() * scale.getr());
    }
    if (++counter == 4)
    {
     for (auto &it : *this)
     {
      it.pseudonormalize();
     }
     counter = 0;
    }
   }
  }
  Base const &normalize() override
  {
   for (auto &it : *this)
   {
    it.normalize();
   }
   return *this;
  }
  base normalize(size_t i) override
  {
   return (*this)[i].normalize();
  }

 private:
  size_t counter = 0;
 };
}


namespace cp_algo::linalg
{
 enum gauss_mode
 {
  normal,
  reverse
 };

 template <typename base_t, class _vec_t = std::conditional_t<
                                math::modint_type<base_t>,
                                modint_vec<base_t>,
                                vec<base_t>>>
 struct matrix : std::vector<_vec_t>
 {
  using vec_t = _vec_t;
  using base = base_t;
  using Base = std::vector<vec_t>;
  using Base::Base;

  matrix(size_t n) : Base(n, vec_t(n)) {}
  matrix(size_t n, size_t m) : Base(n, vec_t(m)) {}

  matrix(Base const &t) : Base(t) {}
  matrix(Base &&t) : Base(std::move(t)) {}

  static matrix from(auto &&r)
  {
   return std::ranges::to<Base>(r);
  }

  size_t n() const { return size(*this); }
  size_t m() const { return n() ? size(row(0)) : 0; }
  auto dim() const { return std::array{n(), m()}; }

  auto &row(size_t i) { return (*this)[i]; }
  auto const &row(size_t i) const { return (*this)[i]; }

  auto operator-() const
  {
   return from(*this | std::views::transform([](auto x)
                                             { return vec_t(-x); }));
  }
  matrix &operator*=(base t)
  {
   for (auto &it : *this)
    it *= t;
   return *this;
  }
  matrix operator*(base t) const { return matrix(*this) *= t; }
  matrix &operator/=(base t) { return *this *= base(1) / t; }
  matrix operator/(base t) const { return matrix(*this) /= t; }

  // Make sure the result is matrix, not Base
  matrix &operator*=(matrix const &t) { return *this = *this * t; }

  void read_transposed()
  {
   for (size_t j = 0; j < m(); j++)
   {
    for (size_t i = 0; i < n(); i++)
    {
     std::cin >> (*this)[i][j];
    }
   }
  }
  void read()
  {
   for (auto &it : *this)
   {
    it.read();
   }
  }
  void print() const
  {
   for (auto const &it : *this)
   {
    it.print();
   }
  }

  static matrix block_diagonal(std::vector<matrix> const &blocks)
  {
   size_t n = 0;
   for (auto &it : blocks)
   {
    assert(it.n() == it.m());
    n += it.n();
   }
   matrix res(n);
   n = 0;
   for (auto &it : blocks)
   {
    for (size_t i = 0; i < it.n(); i++)
    {
     std::ranges::copy(it[i], begin(res[n + i]) + n);
    }
    n += it.n();
   }
   return res;
  }
  static matrix random(size_t n, size_t m)
  {
   matrix res(n, m);
   std::ranges::generate(res, std::bind(vec_t::random, m));
   return res;
  }
  static matrix random(size_t n)
  {
   return random(n, n);
  }
  static matrix eye(size_t n)
  {
   matrix res(n);
   for (size_t i = 0; i < n; i++)
   {
    res[i][i] = 1;
   }
   return res;
  }

  // Concatenate matrices
  matrix operator|(matrix const &b) const
  {
   assert(n() == b.n());
   matrix res(n(), m() + b.m());
   for (size_t i = 0; i < n(); i++)
   {
    res[i] = row(i) | b[i];
   }
   return res;
  }
  matrix submatrix(auto viewx, auto viewy) const
  {
   return from(*this | viewx | std::views::transform([&](auto const &y)
                                                     { return vec_t(y | viewy); }));
  }

  matrix T() const
  {
   matrix res(m(), n());
   for (size_t i = 0; i < n(); i++)
   {
    for (size_t j = 0; j < m(); j++)
    {
     res[j][i] = row(i)[j];
    }
   }
   return res;
  }

  matrix operator*(matrix const &b) const
  {
   assert(m() == b.n());
   matrix res(n(), b.m());
   for (size_t i = 0; i < n(); i++)
   {
    for (size_t j = 0; j < m(); j++)
    {
     res[i].add_scaled(b[j], row(i)[j]);
    }
   }
   return res.normalize();
  }

  vec_t apply(vec_t const &x) const
  {
   return (matrix(1, x) * *this)[0];
  }

  matrix pow(uint64_t k) const
  {
   assert(n() == m());
   return bpow(*this, k, eye(n()));
  }

  matrix &normalize()
  {
   for (auto &it : *this)
   {
    it.normalize();
   }
   return *this;
  }
  template <gauss_mode mode = normal>
  void eliminate(size_t i, size_t k)
  {
   auto kinv = base(1) / row(i).normalize()[k];
   for (size_t j = (mode == normal) * i; j < n(); j++)
   {
    if (j != i)
    {
     row(j).add_scaled(row(i), -row(j).normalize(k) * kinv);
    }
   }
  }
  template <gauss_mode mode = normal>
  void eliminate(size_t i)
  {
   row(i).normalize();
   for (size_t j = (mode == normal) * i; j < n(); j++)
   {
    if (j != i)
    {
     row(j).reduce_by(row(i));
    }
   }
  }
  template <gauss_mode mode = normal>
  matrix &gauss()
  {
   for (size_t i = 0; i < n(); i++)
   {
    eliminate<mode>(i);
   }
   return normalize();
  }
  template <gauss_mode mode = normal>
  auto echelonize(size_t lim)
  {
   return gauss<mode>().sort_classify(lim);
  }
  template <gauss_mode mode = normal>
  auto echelonize()
  {
   return echelonize<mode>(m());
  }

  size_t rank() const
  {
   if (n() > m())
   {
    return T().rank();
   }
   return size(matrix(*this).echelonize()[0]);
  }

  base det() const
  {
   assert(n() == m());
   matrix b = *this;
   b.echelonize();
   base res = 1;
   for (size_t i = 0; i < n(); i++)
   {
    res *= b[i][i];
   }
   return res;
  }

  std::pair<base, matrix> inv() const
  {
   assert(n() == m());
   matrix b = *this | eye(n());
   if (size(b.echelonize<reverse>(n())[0]) < n())
   {
    return {0, {}};
   }
   base det = 1;
   for (size_t i = 0; i < n(); i++)
   {
    det *= b[i][i];
    b[i] *= base(1) / b[i][i];
   }
   return {det, b.submatrix(std::views::take(n()), std::views::drop(n()) | std::views::take(n()))};
  }

  // Can also just run gauss on T() | eye(m)
  // but it would be slower :(
  auto kernel() const
  {
   auto A = *this;
   auto [pivots, free] = A.template echelonize<reverse>();
   matrix sols(size(free), m());
   for (size_t j = 0; j < size(pivots); j++)
   {
    base scale = base(1) / A[j][pivots[j]];
    for (size_t i = 0; i < size(free); i++)
    {
     sols[i][pivots[j]] = A[j][free[i]] * scale;
    }
   }
   for (size_t i = 0; i < size(free); i++)
   {
    sols[i][free[i]] = -1;
   }
   return sols;
  }

  // [solution, basis], transposed
  std::optional<std::array<matrix, 2>> solve(matrix t) const
  {
   matrix sols = (*this | t).kernel();
   if (sols.n() < t.m() || sols.submatrix(
                               std::views::drop(sols.n() - t.m()),
                               std::views::drop(m())) != -eye(t.m()))
   {
    return std::nullopt;
   }
   else
   {
    return std::array{
        sols.submatrix(std::views::drop(sols.n() - t.m()),
                       std::views::take(m())),
        sols.submatrix(std::views::take(sols.n() - t.m()),
                       std::views::take(m()))};
   }
  }

  // To be called after a gaussian elimination run
  // Sorts rows by pivots and classifies
  // variables into pivots and free
  auto sort_classify(size_t lim)
  {
   size_t rk = 0;
   std::vector<size_t> free, pivots;
   for (size_t j = 0; j < lim; j++)
   {
    for (size_t i = rk + 1; i < n() && row(rk)[j] == base(0); i++)
    {
     if (row(i)[j] != base(0))
     {
      std::swap(row(i), row(rk));
      row(rk) = -row(rk);
     }
    }
    if (rk < n() && row(rk)[j] != base(0))
    {
     pivots.push_back(j);
     rk++;
    }
    else
    {
     free.push_back(j);
    }
   }
   return std::array{pivots, free};
  }
 };
 template <typename base_t>
 auto operator*(base_t t, matrix<base_t> const &A) { return A * t; }
}


namespace cp_algo::linalg
{
 enum frobenius_mode
 {
  blocks,
  full
 };
 template <frobenius_mode mode = blocks>
 auto frobenius_form(auto const &A)
 {
  using matrix = std::decay_t<decltype(A)>;
  using vec_t = matrix::vec_t;
  using base = typename matrix::base;
  using base = matrix::base;
  using polyn = math::poly_t<base>;
  assert(A.n() == A.m());
  size_t n = A.n();
  std::vector<polyn> charps;
  std::vector<vec_t> basis, basis_init;
  while (size(basis) < n)
  {
   size_t start = size(basis);
   auto generate_block = [&](auto x)
   {
    while (true)
    {
     vec_t y = x | vec_t::ei(n + 1, size(basis));
     for (auto &it : basis)
     {
      y.reduce_by(it);
     }
     y.normalize();
     if (std::ranges::count(y | std::views::take(n), base(0)) == int(n))
     {
      return polyn(typename polyn::Vector(begin(y) + n, end(y)));
     }
     else
     {
      basis_init.push_back(x);
      basis.push_back(y);
      x = A.apply(x);
     }
    }
   };
   auto full_rec = generate_block(vec_t::random(n));
   // Extra trimming to make it block-diagonal (expensive)
   if constexpr (mode == full)
   {
    if (full_rec.mod_xk(start) != polyn())
    {
     auto charp = full_rec.div_xk(start);
     auto x = basis_init[start];
     auto shift = full_rec / charp;
     for (int j = 0; j < shift.deg(); j++)
     {
      x.add_scaled(basis_init[j], shift[j]);
     }
     basis.resize(start);
     basis_init.resize(start);
     full_rec = generate_block(x.normalize());
    }
   }
   charps.push_back(full_rec.div_xk(start));
  }
  // Find transform matrices while we're at it...
  if constexpr (mode == full)
  {
   for (size_t i = 0; i < n; i++)
   {
    for (size_t j = i + 1; j < n; j++)
    {
     basis[i].reduce_by(basis[j]);
    }
    basis[i].normalize();
   }
   auto T = matrix(basis_init);
   auto Tinv = matrix(basis);
   std::ignore = Tinv.sort_classify(n);
   for (size_t i = 0; i < n; i++)
   {
    Tinv[i] = vec_t(
                  Tinv[i] | std::views::drop(n) | std::views::take(n)) *
              (base(1) / Tinv[i][i]);
   }
   return std::tuple{T, Tinv, charps};
  }
  else
  {
   return charps;
  }
 }

 template <typename base>
 auto with_frobenius(matrix<base> const &A, auto &&callback)
 {
  auto [T, Tinv, charps] = frobenius_form<full>(A);
  std::vector<matrix<base>> blocks;
  for (auto charp : charps)
  {
   matrix<base> block(charp.deg());
   auto xk = callback(charp);
   for (size_t i = 0; i < block.n(); i++)
   {
    std::ranges::copy(xk.a, begin(block[i]));
    xk = xk.mul_xk(1) % charp;
   }
   blocks.push_back(block);
  }
  auto S = matrix<base>::block_diagonal(blocks);
  return Tinv * S * T;
 }

 template <typename base>
 auto frobenius_pow(matrix<base> const &A, uint64_t k)
 {
  return with_frobenius(A, [k](auto const &charp)
                        { return math::poly_t<base>::xk(1).powmod(k, charp); });
 }
};



using namespace std;
using namespace cp_algo::math;

// ---------------------------------------------------------------------------
// Wrapper
// ---------------------------------------------------------------------------
struct CharPoly {
 static vector<u32> run(int n, const vector<vector<u32>>& M_in) {
  using base = cp_algo::math::modint<(int64_t) MOD>;
  using polyn = cp_algo::math::poly_t<base>;
  if (n == 0) return {1};
  cp_algo::linalg::matrix<base> A((size_t) n);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j) A[i][j].setr((uint64_t) M_in[i][j]);
  auto blocks = cp_algo::linalg::frobenius_form(A);
  polyn p = std::reduce(begin(blocks), end(blocks), polyn(1), std::multiplies<polyn>{});
  // p の係数列 (0 次から N 次) を u32 で返す。p_N は 1 (monic)。
  vector<u32> result(n + 1, 0);
  for (int i = 0; i <= n; ++i) {
   if (i < (int) p.a.size()) {
    u32 v = (u32) p.a[i].getr();
    if (v >= MOD) v -= MOD;
    result[i] = v;
   }
  }
  result[n] = 1;
  return result;
 }
};
