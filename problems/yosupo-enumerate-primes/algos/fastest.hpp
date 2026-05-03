#pragma once
#include "_common.hpp"
// yosupo enumerate_primes 最速 (purplesyringa, https://judge.yosupo.jp/submission/336754)
// アルゴリズム部のみ抽出。 I/O は本実験 harness の Solver::run(string) に従い、
// blazingio 等は使用しない。
//
// 全体の流れ (cp_algo::math::sieve_wheel 移植):
//   - 2/3/5/7 wheel (period 210, 48 coprime residues)
//   - bit-packed bit array (u64 word)
//   - dense phase: 小素数群の wheel mask を precompute、 block 内で AND 一発で篩
//   - sparse phase: medium 素数は state machine 風 step で篩
//   - segmented (32MB dense block + 4MB sparse block) で cache 内完結
//   - wheel_primes (= 2,3,5,7) は最終的に sqrt_prime_bits で初期化
//
// big_vector / big_alloc は std::vector に置換 (mmap allocator は不採用)。
// std::ranges::fold_left, std::views::enumerate などは loop で書き直し。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2,bmi,bmi2,popcnt")
#include <immintrin.h>
#endif

#include <cassert>
#include <bit>
#include <new>

namespace yosupo_ep_aligned {
// 32-byte aligned allocator (std::vector default は align 16 程度なので AVX2 で
// vmovdqa が #GP fault する可能性あり、 cp_algo の big_alloc に倣って 32 align)
template<class T>
struct aligned_alloc32 {
 using value_type = T;
 aligned_alloc32() = default;
 template<class U> aligned_alloc32(const aligned_alloc32<U>&) {}
 [[nodiscard]] T* allocate(std::size_t n) {
  return static_cast<T*>(::operator new(n * sizeof(T), std::align_val_t(32)));
 }
 void deallocate(T* p, std::size_t n) noexcept {
  ::operator delete(p, n * sizeof(T), std::align_val_t(32));
 }
 template<class U> bool operator==(const aligned_alloc32<U>&) const { return true; }
 template<class U> bool operator!=(const aligned_alloc32<U>&) const { return false; }
};
}

namespace yosupo_ep_fastest {

using u8 = unsigned char;
using u16 = unsigned short;
using u32 = unsigned int;
using u64 = unsigned long long;

constexpr std::array<u32, 4> wheel_primes = {2u, 3u, 5u, 7u};
constexpr u32 period_v = 2u * 3u * 5u * 7u;       // 210
constexpr u32 coprime_v = 1u * 2u * 4u * 6u;       // φ(210) = 48

constexpr bool coprime_wheel(u32 x) {
 for (u32 p : wheel_primes) if (x % p == 0) return false;
 return true;
}

constexpr auto res_wheel = []() {
 std::array<u8, coprime_v> r{};
 u32 idx = 0;
 for (u32 i = 1; i < period_v; i += 2) if (coprime_wheel(i)) r[idx++] = (u8) i;
 return r;
}();

constexpr auto state_wheel = []() {
 std::array<u8, period_v> s{};
 u8 idx = 0;
 for (u32 i = 0; i < period_v; ++i) {
  s[i] = idx;
  idx += coprime_wheel(i);
 }
 return s;
}();

constexpr auto add_wheel = []() {
 std::array<u8, period_v> a{};
 for (u32 i = 0; i < period_v; ++i) {
  a[i] = 1;
  while (!coprime_wheel(i + a[i])) ++a[i];
 }
 return a;
}();

constexpr auto gap_wheel = []() {
 std::array<u8, coprime_v> g{};
 for (u32 i = 0; i < coprime_v; ++i) g[i] = add_wheel[res_wheel[i]];
 return g;
}();

constexpr u32 to_ord(u32 x) { return (x / period_v) * coprime_v + state_wheel[x % period_v]; }
constexpr u32 to_val(u32 x) { return (x / coprime_v) * period_v + res_wheel[x % coprime_v]; }

// =============================================================================
// bit_array: bit-packed array. static (std::array) and dynamic (std::vector).
// =============================================================================
template<class Cont>
struct _bit_array {
 using word_t = typename Cont::value_type;
 static constexpr size_t width = sizeof(word_t) * 8;
 size_t words = 0, n = 0;
 alignas(32) Cont data{};

 constexpr void resize(size_t N) {
  n = N;
  words = (n + width - 1) / width;
  if constexpr (requires (Cont c) { c.resize(0u); }) data.resize(words);
  else assert(std::size(data) >= words);
 }
 constexpr _bit_array() {
  if constexpr (requires (Cont c) { c.resize(0u); }) {} // dynamic: stays empty until resize
  else resize(std::size(data) * width);
 }
 constexpr _bit_array(size_t N) { resize(N); }
 constexpr word_t& word(size_t x) { return data[x]; }
 constexpr const word_t& word(size_t x) const { return data[x]; }
 constexpr void set_all(word_t v = ~word_t(0)) { for (size_t i = 0; i < words; ++i) data[i] = v; }
 constexpr void reset() { set_all(0); }
 constexpr void set(size_t x) { data[x / width] |= word_t(1) << (x % width); }
 constexpr void reset(size_t x) { data[x / width] &= ~(word_t(1) << (x % width)); }
 constexpr bool test(size_t x) const { return (data[x / width] >> (x % width)) & 1; }
 constexpr bool operator[](size_t x) const { return test(x); }
};

template<size_t N>
using static_bit_array = _bit_array<std::array<u64, (N + 63) / 64>>;
using dynamic_bit_array = _bit_array<std::vector<u64, yosupo_ep_aligned::aligned_alloc32<u64>>>;

template<class BA>
size_t count_bits(const BA& arr, size_t n) {
 size_t res = 0;
 for (size_t i = 0; i < n / BA::width; ++i) res += std::popcount(arr.word(i));
 if (n % BA::width) {
  u64 m = (u64(1) << (n % BA::width)) - 1;
  res += std::popcount(arr.word(n / BA::width) & m);
 }
 return res;
}
template<class BA>
size_t count_bits(const BA& arr) { return count_bits(arr, arr.n); }

// kth set bit in word x (0-indexed)
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
[[gnu::target("bmi2")]]
inline size_t kth_set_bit_in_word(u64 x, size_t k) {
 return std::countr_zero(_pdep_u64(u64(1) << k, x));
}
#else
inline size_t kth_set_bit_in_word(u64 x, size_t k) {
 // fallback: linear scan
 for (size_t pos = 0; pos < 64; ++pos) {
  if ((x >> pos) & 1) {
   if (k == 0) return pos;
   --k;
  }
 }
 return 64;
}
#endif

template<class BA>
size_t skip_bits(const BA& arr, size_t pos, size_t k) {
 if (pos >= arr.n) return arr.n;
 size_t wi = pos / BA::width;
 u64 w = arr.word(wi) >> (pos % BA::width);
 size_t pc = std::popcount(w);
 if (pc > k) return pos + kth_set_bit_in_word(w, k);
 k -= pc;
 while (++wi < arr.words) {
  w = arr.word(wi);
  pc = std::popcount(w);
  if (pc > k) return wi * BA::width + kth_set_bit_in_word(w, k);
  k -= pc;
 }
 return arr.n;
}

// =============================================================================
// sqrt_prime_bits: compile-time bit array of primes ≤ sqrt_threshold
// =============================================================================
constexpr size_t sqrt_threshold = 1 << 15;  // 32768. covers sqrt(N) up to N=2^30
constexpr auto sqrt_prime_bits = []() {
 constexpr int size = sqrt_threshold / 2;
 static_bit_array<size> prime;
 prime.set_all();
 prime.reset(to_ord(1));
 for (u32 i = res_wheel[1]; to_ord(i * i) < size; i += add_wheel[i % period_v]) {
  if (prime[to_ord(i)]) {
   for (u32 k = i; to_ord(i * k) < size; k += add_wheel[k % period_v]) {
    prime.reset(to_ord(i * k));
   }
  }
 }
 return prime;
}();

constexpr size_t num_primes_sqrt = []() {
 size_t cnt = 0;
 for (u32 i = res_wheel[1]; i < sqrt_threshold; i += add_wheel[i % period_v]) {
  cnt += sqrt_prime_bits[to_ord(i)];
 }
 return cnt;
}();
constexpr auto sqrt_primes = []() {
 std::array<u32, num_primes_sqrt> primes{};
 size_t j = 0;
 for (u32 i = res_wheel[1]; i < sqrt_threshold; i += add_wheel[i % period_v]) {
  if (sqrt_prime_bits[to_ord(i)]) primes[j++] = i;
 }
 return primes;
}();

// =============================================================================
// wheel struct + make_wheel
// =============================================================================
struct wheel_t {
 dynamic_bit_array mask;
 u32 product;
};

inline wheel_t make_wheel(const std::vector<u32>& primes, u32 product) {
 assert(product % period_v == 0 && (product / period_v * coprime_v) % dynamic_bit_array::width == 0);
 wheel_t w;
 w.product = product;
 w.mask.resize(product / period_v * coprime_v);
 w.mask.set_all();
 for (u32 p : primes) {
  for (u32 k = 1; p * k < product; k += add_wheel[k % period_v]) {
   w.mask.reset(to_ord(p * k));
  }
 }
 return w;
}

inline void sieve_dense(dynamic_bit_array& prime, u32 l, u32 r, const wheel_t& wheel) {
 if (l >= r) return;
 const u32 width = (u32) dynamic_bit_array::width;
 u32 wl = l / width;
 u32 wr = (r + width - 1) / width;
 u32 N = (u32) wheel.mask.words;
 auto loop = [&](u32 i, u32 block) {
  // assume_aligned<32> は &word(i) で i が 4 の倍数でないと実行時に lie になる
  // (vmovdqa #GP fault)。 unaligned access に任せて vmovdqu を使わせる方が安全。
  u64* p_ptr = &prime.word(i);
  const u64* m_ptr = &wheel.mask.word(0);
  #pragma GCC unroll 48
  for (u32 j = 0; j < block; ++j) p_ptr[j] &= m_ptr[j];
 };
 while (wl + N <= wr) {
  loop(wl, N);
  wl += N;
 }
 if (wl < wr) loop(wl, wr - wl);
}

// sieve_wheel sparse phase: 1 つの medium prime sqrt_primes[i] について
// (l, state) → 篩を打って (new_l, new_state) を返す
inline std::pair<u32, u32> sieve_sparse(dynamic_bit_array& prime, u32 l, u32 r, size_t i, u32 state) {
 // ord_step[i][s] = state s から次の coprime ordinal へのステップ。 size 2*coprime で wrap
 static const auto ord_step = []() {
  std::vector<std::array<u32, 2 * coprime_v>> v(num_primes_sqrt);
  for (u32 i = 0; i < (u32) sqrt_primes.size(); ++i) {
   u32 p = sqrt_primes[i];
   auto& ords = v[i];
   u32 last = to_ord(p);
   for (u32 j = 0; j < coprime_v; ++j) {
    u32 next = to_ord(p * (res_wheel[j] + gap_wheel[j]));
    ords[j] = ords[j + coprime_v] = next - last;
    last = next;
   }
  }
  return v;
 }();
 const auto& ords = ord_step[i];
 u32 p = sqrt_primes[i];
 auto advance = [&]() {
  prime.reset(l);
  l += ords[state++];
 };
 while (l + p * coprime_v <= r) {
  #pragma GCC unroll 48
  for (u32 j = 0; j < coprime_v; ++j) advance();
  state -= coprime_v;
 }
 while (l < r) advance();
 if (state >= coprime_v) state -= coprime_v;
 return {l, state};
}

// メイン: N 以下の素数を bit array で返す
inline dynamic_bit_array sieve_wheel(u32 N) {
 ++N;
 dynamic_bit_array prime(to_ord(N));
 prime.set_all();
 // 各 wheel は小素数群の AND mask。 product を増やしながら作る
 static const auto wheels_data = []() {
  constexpr size_t max_wheel_size = 1 << 20;
  u32 product = period_v * (u32) dynamic_bit_array::width / 4;  // ÷ 2^(|wheel_primes| - 2)
  std::vector<u32> current;
  std::vector<wheel_t> wheels;
  size_t medium_begin = 0;
  for (size_t i = 0; i < sqrt_primes.size(); ++i) {
   u32 p = sqrt_primes[i];
   if (u64(product) * p > max_wheel_size) {
    wheels.push_back(make_wheel(current, product));
    current = {p};
    product = period_v * (u32) dynamic_bit_array::width / 4 * p;
    if (product > max_wheel_size) { medium_begin = i; goto done; }
   } else {
    current.push_back(p);
    product *= p;
   }
  }
  assert(false);
 done:
  return std::pair{wheels, medium_begin};
 }();
 const auto& wheels = wheels_data.first;
 const size_t medium_begin = wheels_data.second;

 // dense phase: 大きい block で AND mask
 constexpr u32 dense_block = 1u << 25;
 for (u32 start = 0; start < N; start += dense_block) {
  u32 r = std::min<u32>(start + dense_block, N);
  for (const auto& wh : wheels) {
   u32 l = start / wh.product * wh.product;
   sieve_dense(prime, to_ord(l), to_ord(r), wh);
  }
 }

 // sparse phase: medium primes のみ別途。 各 prime の (pos, state) を維持
 constexpr u32 sparse_block = 1u << 22;
 std::vector<u32> pos(num_primes_sqrt);
 std::vector<u32> state(num_primes_sqrt);
 for (u32 i = 0; i < (u32) sqrt_primes.size(); ++i) {
  u32 p = sqrt_primes[i];
  pos[i] = to_ord(p * p);
  state[i] = state_wheel[p % period_v];
 }
 for (u32 start = 0; start < N; start += sparse_block) {
  u32 r = to_ord(std::min<u32>(start + sparse_block, N));
  for (size_t i = medium_begin; i < sqrt_primes.size(); ++i) {
   if (state[i] >= r) break;  // pos[i] が r を超えたら stop... ※ original logic では state[i] というより pos[i] と r の比較が自然だが移植元準拠
   if (pos[i] >= r) break;
   auto [np, ns] = sieve_sparse(prime, pos[i], r, i, state[i]);
   pos[i] = np;
   state[i] = ns;
  }
 }
 // sqrt_prime_bits を最初の word 群にコピー (wheel_primes も含めて prime 判定を正す)
 for (size_t i = 0; i < std::min(prime.words, sqrt_prime_bits.words); ++i) {
  prime.word(i) = sqrt_prime_bits.word(i);
 }
 return prime;
}

}  // namespace

struct Solver {
 static std::pair<u32, std::vector<u32>> run(u32 N, u32 A, u32 B) {
  using namespace yosupo_ep_fastest;

  auto primes_bits = sieve_wheel(N);
  size_t cnt = count_bits(primes_bits);
  size_t extra = 0;
  for (u32 p : wheel_primes) if (p <= N) ++extra;
  cnt += extra;
  size_t X = (cnt < B) ? 0 : (cnt - B + A - 1) / A;
  std::vector<u32> selected;
  selected.reserve(X);

  size_t b_remaining = B;
  size_t x_remaining = X;
  for (u32 p : wheel_primes) {
   if (b_remaining == 0 && x_remaining > 0 && p <= N) {
    selected.push_back(p);
    --x_remaining;
   }
   if (p <= N) {
    if (b_remaining > 0) --b_remaining;
    else b_remaining = A - 1;
   }
  }
  size_t pos = skip_bits(primes_bits, 0, b_remaining);
  while (pos < primes_bits.n && x_remaining > 0) {
   u32 v = to_val((u32) pos);
   if (v > N) break;
   selected.push_back(v);
   --x_remaining;
   if (x_remaining > 0) pos = skip_bits(primes_bits, pos + 1, A - 1);
  }
  return {(u32) cnt, std::move(selected)};
 }
};
