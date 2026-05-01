#pragma once
// =============================================================================
// 参考: yosupo "Primality Test" 提出 139990 (Rust, mizar の cargo 由来).
//   https://judge.yosupo.jp/submission/139990
//   方式: Baillie-PSW (BPSW) primality test。
//     is_prime(n) := miller_rabin(n, base=2) && strong_lucas(n)
//   テーブル lookup は使わず、Lucas 数列 + Jacobi 記号を 1 回計算するだけ。
//   現時点で 64-bit 全域に counter-example は知られていない (実用上決定的)。
//   Lucas test 用のパラメータ D は Selfridge 法: D = 5, -7, 9, -11, 13, ...
//   と試して Jacobi(D, n) = -1 になる最初の D を採用。
// 抽出方針:
//   - 元は Rust の cargo bundle なのでロジックを C++ に port
//   - Mont (Montgomery64), Jacobi, Lucas, MR(2) を最小限実装
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"

namespace yosupo_bpsw {
using std::uint32_t;
using std::uint64_t;
using std::int64_t;
using u128_local = __uint128_t;

// Montgomery 64. n は奇数前提。R = 2^64。
// 慣習: mi は n * mi ≡ +1 (mod 2^64) を満たす方 (Rust 元実装と同じ向き)。
struct Mont {
 uint64_t m;     // modulus (odd)
 uint64_t mi;    // n * mi ≡ 1 mod 2^64
 uint64_t mh;    // (m+1)/2  (used for div2)
 uint64_t R;     // R mod m
 uint64_t Rn;    // -R mod m = m - R
 uint64_t R2;    // R^2 mod m
 uint64_t d;     // (m-1) >> k
 int k;          // ctz(m-1)
 explicit Mont(uint64_t mm): m(mm) {
  // Newton iteration for mi: m * mi ≡ 1 (mod 2^64).
  uint64_t x = mm * 3 ^ 2;
  uint64_t y = (uint64_t) 1 - mm * x;
  for (int w = 5; w < 64; w *= 2) { x = x * (y + 1); y = y * y; }
  mi = x;
  mh = (mm >> 1) + 1;
  R = (uint64_t) -(int64_t) 1 % mm + 1;     // 2^64 mod m
  R2 = (uint64_t) ((u128_local) -(__int128) 1 % mm) + 1; // 2^128 mod m
  Rn = mm - R;
  d = mm - 1;
  k = __builtin_ctzll(d);
  d >>= k;
 }
 inline uint64_t add_mod(uint64_t a, uint64_t b) const {
  uint64_t t = a - (m - b);
  return t + ((-(int64_t)(b > m || a + b < a) | (a + b >= m)) ? 0 : 0); // 単純化
 }
 inline uint64_t add(uint64_t a, uint64_t b) const {
  // a, b < m。a + b mod m。
  uint64_t s = a + b;
  return s >= m ? s - m : s;
 }
 inline uint64_t sub(uint64_t a, uint64_t b) const {
  return a >= b ? a - b : a + (m - b);
 }
 inline uint64_t div2(uint64_t ar) const {
  // (ar / 2) mod m。ar が偶数ならそのまま、奇数なら +m してから 2 で割る。
  return (ar >> 1) + ((ar & 1) ? mh : 0);
 }
 // Montgomery multiply: (ar * br * R^{-1}) mod m。
 inline uint64_t mrmul(uint64_t ar, uint64_t br) const {
  u128_local t = (u128_local) ar * br;
  uint64_t q = (uint64_t) t * mi;
  uint64_t hi = (uint64_t) (t >> 64);
  uint64_t mq_hi = (uint64_t) (((u128_local) q * m) >> 64);
  uint64_t res = hi - mq_hi;
  return hi < mq_hi ? res + m : res;
 }
 // 通常値 a → Montgomery 形式 ar = a*R mod m。
 inline uint64_t to_mont(uint64_t a) const { return mrmul(a, R2); }
 // Mont 形式 ar → 通常値 a = ar * R^{-1} mod m
 inline uint64_t from_mont(uint64_t ar) const { return mrmul(ar, 1); }
 inline uint64_t mod_n(uint64_t x) const { return x < m ? x : x % m; }
 inline uint64_t pow_mont(uint64_t ar, uint64_t b) const {
  if (b == 0) return R;
  while ((b & 1) == 0) { ar = mrmul(ar, ar); b >>= 1; }
  uint64_t t = ar; b >>= 1;
  while (b) {
   ar = mrmul(ar, ar);
   if (b & 1) t = mrmul(t, ar);
   b >>= 1;
  }
  return t;
 }
};

// Miller-Rabin with base 2 (strong probable prime test).
inline bool mr_base2(const Mont& mo) {
 // base 2 in Mont = R*2 mod m = add(R, R)
 uint64_t br = mo.add(mo.R, mo.R);
 br = mo.pow_mont(br, mo.d);
 if (br == mo.R || br == mo.Rn) return true;
 for (int i = 1; i < mo.k; ++i) {
  br = mo.mrmul(br, br);
  if (br == mo.Rn) return true;
 }
 return false;
}

// Jacobi 記号 (a/n) for u64 a (positive) and odd n.
// 戻り値: -1, 0, 1。
inline int jacobi_u64(uint64_t a, uint64_t n) {
 int j = 1;
 while (a) {
  int b = __builtin_ctzll(a);
  if ((b & 1) && ((n + 2) & 5) == 5) j = -j;
  a >>= b;
  if (a < n) {
   if (a > (n >> 1)) {
    a = n - a;
    if ((n & 3) == 3) j = -j;
    int b2 = __builtin_ctzll(a);
    if ((b2 & 1) && ((n + 2) & 5) == 5) j = -j;
    a >>= b2;
   }
   if ((a & n & 3) == 3) j = -j;
   uint64_t t = n - a;
   n = a;
   a = t;
  } else {
   a -= n;
  }
 }
 return n == 1 ? j : 0;
}
// 整数値 (符号付き) の Jacobi。a < 0 のとき (-1/n) を補正。
inline int jacobi_signed(int64_t ai, uint64_t n) {
 if (ai >= 0) return jacobi_u64((uint64_t) ai, n);
 int sign = ((n & 3) == 3) ? -1 : 1;
 return sign * jacobi_u64((uint64_t) -ai, n);
}

// Selfridge-style D selection: D = 5, -7, 9, -11, 13, ...
//   d = 5; while jacobi(d, n) != -1: d = -(d+2*sign(d))*sign(d) ... ?
// Rust の lucas_dprobe を直訳:
//   for d in (5u32..).step_by(2):
//     j = jacobi((n%4==3 && d%4==3) ? -jacobi(d,n) : jacobi(d,n))
//     match j: -1 -> return d (or -d if d&2==2)
//              0  -> if d < n return None  (= n is composite, d divides n)
//              1  -> continue
//   plus: 試行回数 30 程度を超えたら平方数チェックで打ち切り。
// 簡略化: d=5,-7,9,-11,... を直接列挙。
inline bool is_perfect_square(uint64_t n) {
 if (n == 0) return true;
 uint64_t s = (uint64_t) std::sqrt((double) n);
 for (uint64_t t = (s ? s - 1 : 0); t <= s + 1; ++t) if (t * t == n) return true;
 return false;
}
inline int lucas_dprobe(uint64_t n) {
 // d = 5, -7, 9, -11, 13, -15, ...  symbolically; signed.
 for (uint32_t mag = 5;; mag += 2) {
  int sign = (mag & 2) ? -1 : 1; // mag=5 → 5, 7→-7, 9→9, 11→-11
  int64_t D = (int64_t) mag * sign;
  int j = jacobi_signed(D, n);
  if (j == -1) return (int) D;
  if (j == 0 && mag < n) return 0; // n は D の倍数 → 合成数
  // ロックイン: 何度試しても -1 にならないなら平方数チェックで早期終了。
  if (mag == 61) {
   if (is_perfect_square(n)) return 0;
  }
  if (mag > 200000) return 0; // safety net
 }
}

// Strong Lucas pseudoprime test (P=1, Q=(1-D)/4)。
// 元 Rust 実装を直訳。
inline bool lucas_test(const Mont& mo) {
 uint64_t n = mo.m;
 int D = lucas_dprobe(n);
 if (D == 0) return false;
 // Q = (1 - D) / 4 if D < 0 else n - (D - 1)/4 (= -(D-1)/4 mod n)
 uint64_t Q;
 if (D < 0) Q = (uint64_t) (1 - D) / 4;
 else Q = n - (uint64_t) (D - 1) / 4;
 uint64_t qm = mo.to_mont(mo.mod_n(Q));
 uint64_t dm = (D < 0) ? mo.to_mont(n - mo.mod_n((uint64_t) -D))
                       : mo.to_mont(mo.mod_n((uint64_t) D));

 uint64_t um = mo.R; // U_1 = 1
 uint64_t vm = mo.R; // V_1 = P = 1
 uint64_t qn = qm;   // Q^1
 // 主ループ: bits of (n+1) を MSB から処理 (leading 1 はスキップ済み)。
 uint64_t k = (n + 1) << __builtin_clzll(n + 1);
 k <<= 1; // skip leading 1
 while (k) {
  // doubling: U_{2m} = U_m V_m, V_{2m} = V_m^2 - 2 Q^m
  um = mo.mrmul(um, vm);
  vm = mo.sub(mo.mrmul(vm, vm), mo.add(qn, qn));
  qn = mo.mrmul(qn, qn);
  if ((int64_t) k < 0) { // top bit set → 1-step
   // U_{m+1} = (U_m + V_m) / 2, V_{m+1} = (D U_m + V_m) / 2  (with P=1)
   uint64_t uu = mo.div2(mo.add(um, vm));
   vm = mo.div2(mo.add(mo.mrmul(dm, um), vm));
   um = uu;
   qn = mo.mrmul(qn, qm);
  }
  k <<= 1;
 }
 // 主ループ終了時: (um, vm, qn) = (U_d, V_d, Q^d)。
 if (um == 0 || vm == 0) return true;
 // Strong test: V_{d * 2^r} == 0 を r=1..s でチェック。
 uint64_t x = (n + 1) & ~n; // = 2^s where s = ctz(n+1)
 x >>= 1;
 while (x) {
  um = mo.mrmul(um, vm);
  vm = mo.sub(mo.mrmul(vm, vm), mo.add(qn, qn));
  if (vm == 0) return true;
  qn = mo.mrmul(qn, qn);
  x >>= 1;
 }
 return false;
}

inline bool is_prime(uint64_t n) {
 if (n < 4) return n >= 2;
 if ((n & 1) == 0) return false;
 // 小さい素数で trial division (高速化)
 for (uint64_t p : {3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull, 41ull, 43ull})
  if (n % p == 0) return n == p;
 Mont mo(n);
 return mr_base2(mo) && lucas_test(mo);
}
} // namespace yosupo_bpsw

struct Primality {
 static vector<bool> run(const vector<u64>& qs) {
  vector<bool> ans(qs.size());
  for (size_t i = 0; i < qs.size(); ++i)
   ans[i] = yosupo_bpsw::is_prime(qs[i]);
  return ans;
 }
};
