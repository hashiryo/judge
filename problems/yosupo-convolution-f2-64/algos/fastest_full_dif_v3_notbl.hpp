#pragma once
#include "_common.hpp"
// fastest_full_dif_v3.hpp の **submask table 不要版**。
// (v3 = v2 + 4-way unroll、 ここでは加えて SUBMASK_TBL を撤去。)
//
// アルゴリズム:
//   - Phase A (true DIF, no shuffle): monomial → LCH 自然順
//     再帰的 P = Q_0 + s_{k-1}(x) · Q_1 分解 (Cantor subspace poly)。
//     Cantor で s_k(x) = Σ_{l ⊆ k bits} x^{2^l}。
//   - Phase B (DIF top-down butterflies, single ska per block):
//     LCH natural → eval natural。 ska は level 非依存な master 表 M[j] から。
//
// DIT 全部入り (fastest_full_v2) との対応:
//
// 1. Cantor 対称性 (DIT 用 2 mul → 1 mul):
//    DIF butterfly は元から (lo + ska·hi, lo + (ska+1)·hi) = (t, t ^ hi) で 1 mul。
//    DIT で言う Cantor 対称性を構造的に吸収済み。 移植不要。
//
// 2. halftable (twiddle メモリ半分):
//    DIF master M[j] = Σ_{L: j_L=1} β_{L+1} は元から size 2^(d-1)、 level 非依存
//    で 1 本だけ。 DIT 全部入りの per-level halftable よりさらにコンパクト。
//
// 3. ska=0 ブロック完全 skip + mul2 ペアリング:
//    block j=0 は ska=0 → ブロック全体が butterfly_0 (PCLMUL 無し) で済む。
//    j>0 のみ btf_fft_block_dif で SIMD 処理。
//
// 4. butterfly pair の SIMD pipeline (btf_fft_block_dif):
//    block 内の (i, i+1) ペアを __m256i で。 ska は scalar、 一度 broadcast して
//    block 全体で再利用 (DIT は g0[i] が i 毎に変わるのでロードが必要だったが
//    DIF はロード 1 回のみ)。 それだけ register pressure が低く pipeline がきれい。
//
// 5. Phase A 最終段 fusion:
//    DIF Phase A の総 op 数は DIT m-loop の ~1/5 (Σpopcount(k-1) vs Σ(k-1))。
//    levels 2, 3 (len=4, 8) は inner ops が単純 (popcount(1)=1, popcount(2)=1)、
//    特殊化して unroll する。
//
// 6. init() 簡略化: 既達。 master 1 本だけ。
//
// 必要拡張: PCLMUL + VPCLMULQDQ + AVX2 + BMI2 (Ice Lake / Zen3 以降)。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
#include "../../gf2-64/_shared/sq.hpp"

namespace conv_f2_64_full_dif_v3_notbl {

using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// =============================================================================
// constexpr GF(2^64) mul / sq
// =============================================================================
constexpr void clmul128_ce(u64 a, u64 b, u64& lo_out, u64& hi_out) {
 u64 Tlo[16]= {0}, Thi[16]= {0};
 Tlo[1]= a;
 for(int v= 2; v < 16; ++v) {
  u64 plo= Tlo[v >> 1], phi= Thi[v >> 1];
  u64 nlo= plo << 1;
  u64 nhi= (phi << 1) | (plo >> 63);
  if(v & 1) nlo^= a;
  Tlo[v]= nlo;
  Thi[v]= nhi;
 }
 u64 lo= 0, hi= 0;
 for(int s= 60; s >= 0; s-= 4) {
  u64 nhi= (hi << 4) | (lo >> 60);
  u64 nlo= lo << 4;
  u32 nib= u32((b >> s) & 0xF);
  lo= nlo ^ Tlo[nib];
  hi= nhi ^ Thi[nib];
 }
 lo_out= lo;
 hi_out= hi;
}
constexpr u64 mul_ce(u64 a, u64 b) {
 u64 lo, hi;
 clmul128_ce(a, b, lo, hi);
 u64 f1l, f1h, f2l, f2h;
 clmul128_ce(hi, 0x1B, f1l, f1h);
 clmul128_ce(f1h, 0x1B, f2l, f2h);
 return lo ^ f1l ^ f2l;
}
constexpr u64 sq_ce(u64 a) { return mul_ce(a, a); }

// =============================================================================
// Cantor chain。 β_l = CHAIN[62-l]、 β_l^2 + β_l = β_{l-1}、 β_0 = 1
// =============================================================================
constexpr int CHAIN_LEN= 63;
constexpr auto CHAIN= []() {
 array<u64, CHAIN_LEN> c{};
 c[0]= 2;
 for(int k= 1; k < CHAIN_LEN; ++k) c[k]= sq_ce(c[k - 1]) ^ c[k - 1];
 return c;
}();

template<class T> int msb(T n) { return n == 0 ? -1 : 63 - __builtin_clzll(n); }
template<class T> T ceil_pow2(T n) { return n <= 1 ? T(1) : T(1) << (msb(n - 1) + 1); }

// =============================================================================
// DIF master twiddle (single contiguous, level 非依存)
//   M[j] = Σ_{L: j_L=1} β_{L+1}   (j ∈ [0, 2^(d-1)))
// =============================================================================
struct nim_fft_data {
 std::vector<u64> master;
 int max_d= -1;
 void init(int d) {
  if(d <= max_d) return;
  max_d= d;
  int sz= (d <= 1) ? 1 : (1 << (d - 1));
  master.assign(sz, 0);
  for(int j= 1; j < sz; ++j) {
   int L= __builtin_ctz(j);
   master[j]= master[j ^ (1 << L)] ^ CHAIN[61 - L];
  }
 }
};
inline nim_fft_data nim_data;

// =============================================================================
// SIMD primitives (fastest_full_v2 から借用、 ska broadcast 化)
// =============================================================================
VPCLMUL inline __m256i expand_pair(const u64* p) {
 __m128i v= _mm_loadu_si128((const __m128i*)p);
 return _mm256_permute4x64_epi64(_mm256_castsi128_si256(v), _MM_SHUFFLE(3, 1, 2, 0));
}
VPCLMUL inline __m128i pack_pair(__m256i v) {
 return _mm256_castsi256_si128(_mm256_permute4x64_epi64(v, _MM_SHUFFLE(3, 2, 2, 0)));
}
inline const __m256i RED_TABLE_DIF3= _mm256_setr_epi8(
 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);

VPCLMUL inline __m256i clmul_reduce_pair(__m256i a_vec, __m256i b_vec) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i h_idx= _mm256_srli_epi64(prod, 60);
 __m256i indices= _mm256_srli_si256(h_idx, 8);
 __m256i red_vec= _mm256_shuffle_epi8(RED_TABLE_DIF3, indices);
 return _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), red_vec);
}

// DIF FFT butterfly (block 内 ska 一定): pair (i, i+1) を SIMD で。
//   r = mul(hi, ska); lo' = lo ^ r; hi' = lo' ^ hi
// ska_vec は呼び出し側で 1 度だけ broadcast されたもの。
VPCLMUL inline void btf_fft_pair_dif(u64* f0_ptr, u64* f1_ptr, __m256i ska_vec) {
 __m256i hi_vec= expand_pair(f1_ptr);
 __m256i p_vec= clmul_reduce_pair(hi_vec, ska_vec);
 __m256i lo_vec= expand_pair(f0_ptr);
 __m256i t_vec= _mm256_xor_si256(lo_vec, p_vec);
 __m256i f1_new= _mm256_xor_si256(t_vec, hi_vec);
 _mm_storeu_si128((__m128i*)f0_ptr, pack_pair(t_vec));
 _mm_storeu_si128((__m128i*)f1_ptr, pack_pair(f1_new));
}

// DIF IFFT butterfly: 順序逆。
//   hi_new = lo ^ hi (= s); lo_new = lo ^ s*ska
VPCLMUL inline void btf_ifft_pair_dif(u64* f0_ptr, u64* f1_ptr, __m256i ska_vec) {
 __m256i lo_vec= expand_pair(f0_ptr);
 __m256i hi_vec= expand_pair(f1_ptr);
 __m256i s_vec= _mm256_xor_si256(lo_vec, hi_vec);
 __m256i p_vec= clmul_reduce_pair(s_vec, ska_vec);
 __m256i lo_new= _mm256_xor_si256(lo_vec, p_vec);
 _mm_storeu_si128((__m128i*)f0_ptr, pack_pair(lo_new));
 _mm_storeu_si128((__m128i*)f1_ptr, pack_pair(s_vec));
}

// ska scalar → __m256i 形式 (lane 0, 2 のみ valid、 pair clmul の入力レイアウト)
// 構造: (ska, 0, ska, 0)
VPCLMUL inline __m256i broadcast_ska(u64 ska) {
 __m128i lo= _mm_set_epi64x(0, ska);
 __m256i v= _mm256_castsi128_si256(lo);
 return _mm256_inserti128_si256(v, lo, 1);
}

// =============================================================================
// Phase A (true DIF): bc_to_lch / bc_to_mono
//   level 2, 3 は popcount(level-1) = 1 で単純 (1 submask {0})、 inner ops 1/i。
//   level k ≥ 4 は popcount(k-1) ≥ 1 (一般に複数 submask)。
// =============================================================================
inline void bc_to_lch(u64* poly, int n) {
 if(n <= 1) return;
 int d= msb(n);
 for(int level= d; level >= 2; --level) {
  int len= 1 << level;
  int half= len >> 1;
  int sub_idx= level - 1;
  for(int base= 0; base < n; base+= len) {
   u64* p= poly + base;
   for(int i= half - 1; i >= 0; --i) {
    u64 q= p[half + i];
    if(q == 0) continue;
    int s= sub_idx;
    while(true) {
     s= (s - 1) & sub_idx;
     p[i + (1 << s)]^= q;
     if(s == 0) break;
    }
   }
  }
 }
}

inline void bc_to_mono(u64* poly, int n) {
 if(n <= 1) return;
 int d= msb(n);
 for(int level= 2; level <= d; ++level) {
  int len= 1 << level;
  int half= len >> 1;
  int sub_idx= level - 1;
  for(int base= 0; base < n; base+= len) {
   u64* p= poly + base;
   for(int i= 0; i < half; ++i) {
    u64 q= p[half + i];
    if(q == 0) continue;
    int s= sub_idx;
    while(true) {
     s= (s - 1) & sub_idx;
     p[i + (1 << s)]^= q;
     if(s == 0) break;
    }
   }
  }
 }
}

// =============================================================================
// Phase B: DIF top-down butterflies、 j=0 完全 skip + i pair SIMD
// =============================================================================
VPCLMUL inline void nim_fft(std::vector<u64>& f) {
 int n= (int)f.size();
 if(n <= 1) return;
 int d= msb(n);
 bc_to_lch(f.data(), n);
 const u64* M= nim_data.master.data();
 for(int k= d; k >= 1; --k) {
  int unit= 1 << k;
  int half= unit >> 1;
  int num_blocks= n / unit;
  // j=0: ska=0 で butterfly_0 (PCLMUL 不要)
  {
   u64* base= f.data();
   for(int i= 0; i < half; ++i) base[half + i]^= base[i];
  }
  // j>0: ska scalar broadcast、 i pair SIMD
  for(int j= 1; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * unit;
   u64 ska= M[j];
   if(half == 1) {
    // SIMD pair が立てられない (1 個しかない)
    u64 r= mul(base[1], ska);
    base[0]^= r;
    base[1]^= base[0];
    continue;
   }
   __m256i ska_vec= broadcast_ska(ska);
   int i= 0;
   // 4 個の VPCLMULQDQ を独立に発行 → OoO で latency hiding
   for(; i + 7 < half; i+= 8) {
    btf_fft_pair_dif(base + i, base + half + i, ska_vec);
    btf_fft_pair_dif(base + i + 2, base + half + i + 2, ska_vec);
    btf_fft_pair_dif(base + i + 4, base + half + i + 4, ska_vec);
    btf_fft_pair_dif(base + i + 6, base + half + i + 6, ska_vec);
   }
   for(; i + 1 < half; i+= 2) {
    btf_fft_pair_dif(base + i, base + half + i, ska_vec);
   }
   for(; i < half; ++i) {
    u64 r= mul(base[half + i], ska);
    base[i]^= r;
    base[half + i]^= base[i];
   }
  }
 }
}

VPCLMUL inline void nim_ifft(std::vector<u64>& f) {
 int n= (int)f.size();
 if(n <= 1) return;
 int d= msb(n);
 const u64* M= nim_data.master.data();
 for(int k= 1; k <= d; ++k) {
  int unit= 1 << k;
  int half= unit >> 1;
  int num_blocks= n / unit;
  // j=0
  {
   u64* base= f.data();
   for(int i= 0; i < half; ++i) base[half + i]^= base[i];
  }
  for(int j= 1; j < num_blocks; ++j) {
   u64* base= f.data() + (size_t)j * unit;
   u64 ska= M[j];
   if(half == 1) {
    base[1]^= base[0];
    base[0]^= mul(base[1], ska);
    continue;
   }
   __m256i ska_vec= broadcast_ska(ska);
   int i= 0;
   for(; i + 7 < half; i+= 8) {
    btf_ifft_pair_dif(base + i, base + half + i, ska_vec);
    btf_ifft_pair_dif(base + i + 2, base + half + i + 2, ska_vec);
    btf_ifft_pair_dif(base + i + 4, base + half + i + 4, ska_vec);
    btf_ifft_pair_dif(base + i + 6, base + half + i + 6, ska_vec);
   }
   for(; i + 1 < half; i+= 2) {
    btf_ifft_pair_dif(base + i, base + half + i, ska_vec);
   }
   for(; i < half; ++i) {
    base[half + i]^= base[i];
    base[i]^= mul(base[half + i], ska);
   }
  }
 }
 bc_to_mono(f.data(), n);
}

VPCLMUL inline std::vector<u64> nim_convolution(std::vector<u64> f, std::vector<u64> g) {
 int n= (int)f.size(), m= (int)g.size();
 int s= (int)ceil_pow2(u32(n + m - 1));
 f.resize(s);
 g.resize(s);
 nim_data.init(msb(s));
 nim_fft(f);
 nim_fft(g);
 for(int i= 0; i < s; ++i) f[i]= mul(f[i], g[i]);
 nim_ifft(f);
 f.resize(n + m - 1);
 return f;
}

}  // namespace conv_f2_64_full_dif_v3_notbl

struct Solver {
 VPCLMUL static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_full_dif_v3_notbl;
  auto c= nim_convolution(a_in, b_in);
  return c;
 }
};
