#pragma once
// pclmul_pohlig_v11 + Packed BSGS + Perfect-hash 641 + Packed 65537:
//
// 改良点 (v11 比):
//   1. BSGSTable6700417: SoA (keys+values) -> Packed u64 entries (B1 layout).
//      - 1 lookup = 1 cache line (was 2).
//      - inv_base_m を pre-computed constexpr 定数に置換 (build 時の pow_bw 撤廃).
//      - if(res < q) ガードは構造的に res >= q が起きないため削除.
//   2. direct_65537: SoA -> Packed u32 entries (15-bit fingerprint + 17-bit value).
//      - cap 倍増 (262144 -> 524288, α=0.25 -> 0.125) で probe 削減.
//      - メモリ 3 MB -> 2 MB.
//   3. direct_641: 完全ハッシュ (Fibonacci-style multiplicative hash).
//      - 衝突ゼロなので probe / EMPTY check 不要.
//      - constexpr build (デフォルト constexpr limits 内).
//      - .rodata 焼き込み, init_tables() から削除.
//      - メモリ 24 KB -> 32 KB, lookup 3.88 ns -> 0.55 ns.
//
// 必要な拡張: VPCLMULQDQ + AVX2 (Intel Ice Lake / AMD Zen3 以降).
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_log_pohlig_v11 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob2;
using gf2_64_pclmul::frob3;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob6;
using gf2_64_pclmul::frob7;
using gf2_64_pclmul::frob8;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
// VPCLMUL 2 並列 mul + 並列 reduction (vmul_3_2 と同じ idiom)
VPCLMUL inline __m256i mul2(__m256i a_vec, __m256i b_vec, u64& r0, u64& r1) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i h_idx= _mm256_srli_epi64(prod, 60);
 __m256i indices= _mm256_srli_si256(h_idx, 8);
 __m256i red_vec= _mm256_shuffle_epi8(RED_TABLE, indices);
 __m256i result= _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), red_vec);
 r0= _mm256_extract_epi64(result, 0);
 r1= _mm256_extract_epi64(result, 2);
 return result;
}
// =============================================================================
// constexpr GF(2^64) 乗算 (PerfectHash641 の build 用).
// PCLMUL intrinsic は constexpr 化できないため、4-bit windowed CLMUL で実装。
// reduction polynomial は runtime 版と同一: x^64 + x^4 + x^3 + x + 1 (R = 0x1B).
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
 constexpr u64 R= 0x1B;
 u64 fold1_lo, fold1_hi;
 clmul128_ce(hi, R, fold1_lo, fold1_hi);
 u64 fold2_lo, fold2_hi;
 clmul128_ce(fold1_hi, R, fold2_lo, fold2_hi);
 return lo ^ fold1_lo ^ fold2_lo;
}
// =============================================================================
// F_{2^16}^* log table (compile-time 構築) - 元の v11 から変更なし
// =============================================================================
constexpr auto LN16= []() {
 u16 col[]= {1U, 11778U, 7028U, 51115U, 48663U, 26081U, 17458U, 40223U, 30334U, 42368U, 14380U, 2223U, 49688U, 11217U, 44239U, 63445U};
 u16 T_lo[256]= {}, T_hi[256]= {};
 for(int v= 0; v < 256; ++v) {
  u16 lo= 0, hi= 0;
  for(int j= 0; j < 8; ++j)
   if((v >> j) & 1) {
    lo^= col[j];
    hi^= col[j + 8];
   }
  T_lo[v]= lo;
  T_hi[v]= hi;
 }
 array<u16, 65536> ln{};
 u16 cur= 1;
 for(u32 k= 0; k < 65535; ++k) {
  u16 lo= T_lo[u8(cur)] ^ T_hi[cur >> 8];
  ln[lo]= u32(u16(k)) * 2699 % 65535;
  cur= u16(cur << 1) ^ (0x002DU & -u16(cur >> 15));
 }
 ln[0]= 0;
 return ln;
}();
// =============================================================================
// 定数群
// =============================================================================
constexpr u64 P_F16= 65535;
constexpr u64 P_641= 641;
constexpr u64 P_F17= 65537;
constexpr u64 P_BIG= 6700417;
constexpr u64 INV_65535_641= 243ULL;
constexpr u64 MOD2= (P_F16 * P_641) % P_F17;
constexpr u64 INV_MOD2_F17= 45242ULL;
constexpr u64 MOD3= (P_F16 * P_641 * P_F17) % P_BIG;
constexpr u64 INV_MOD3_BIG= 3883315ULL;
constexpr u64 MOD_F16= P_F16;
constexpr u64 MOD_F16_641= P_F16 * P_641;
constexpr u64 MOD_F16_641_F17= P_F16 * P_641 * P_F17;
constexpr u64 EXP_F16= 0x0001000100010001ull;
constexpr u64 EXP_641= 0x00663d80ff99c27full;
constexpr u64 EXP_F17= 0x0000ffff0000ffffull;
constexpr u64 EXP_BIG= 0x00000280fffffd7full;
constexpr u64 G_641= 0x6bf808f7824282a2ull;
constexpr u64 G_65537= 0x1c1e79669b95a7ceull;
constexpr u64 G_6700417= 0x00f542601703f991ull;
// =============================================================================
// PerfectHash641: 641 元 subgroup 用 完全ハッシュテーブル (constexpr).
//
// ハッシュ関数: h(x) = (x * C) >> (64 - BITS), Fibonacci-style.
// C = 0xffef5fb99f1bf6e7 は ~50 万ランダム試行で衝突ゼロを発見した定数。
// cap = 2^14 = 16384, u16 entries -> table size = 32 KB.
//
// lookup target は必ず order-641 subgroup の元なので、未登録 key を引くことはない。
// したがって empty marker / fp check は不要、tab[hash(t)] を直接返すだけ。
// =============================================================================
struct PerfectHash641 {
 static constexpr u64 C= 0xffef5fb99f1bf6e7ull;
 static constexpr u32 BITS= 14;
 static constexpr u32 CAP= 1u << BITS;
 static constexpr u32 SHIFT= 64 - BITS;
 static constexpr array<u16, CAP> tab= []() {
  array<u16, CAP> t{};
  u64 cur= 1;
  for(u32 k= 0; k < P_641; ++k) {
   t[u32((cur * C) >> SHIFT)]= u16(k);
   cur= mul_ce(cur, G_641);
  }
  return t;
 }();
 static u32 lookup(u64 t) { return tab[u32((t * C) >> SHIFT)]; }
};
// =============================================================================
// DirectLogTable_65537: 65537 元 subgroup 用 packed linear-probing.
//
// レイアウト: u32 entry = (fp << 17) | (k + 1)
//   - 下位 17 bit: k + 1 (k ∈ [0, 65537), empty = 0 と区別するため +1)
//   - 上位 15 bit: key の fp = (key >> BITS) & ((1 << 15) - 1)
//
// fp 15 bit + k+1 17 bit = 32 bit にぴったり収まる。
// cap = 2^19 = 524288 (load factor α = 0.125), メモリ = 2 MB.
//
// 全数 lookup チェックで衝突ゼロを確認済み (生成元固定なのでテーブルも固定).
// =============================================================================
struct DirectLogTable_65537 {
 static constexpr u32 BITS= 19;
 static constexpr u32 MASK= (1u << BITS) - 1;
 static constexpr u32 CAP= MASK + 1;
 static constexpr u32 VBITS= 17;
 static constexpr u32 VMASK= (1u << VBITS) - 1;
 static constexpr u32 FP_MASK= (1u << (32 - VBITS)) - 1;  // 15 bit
 static constexpr u32 EMPTY= 0;
 std::vector<u32> tab;
 void build() {
  tab.assign(CAP, EMPTY);
  u64 cur= 1;
  for(u32 k= 0; k < P_F17; ++k) {
   u32 h= u32(cur) & MASK;
   while(tab[h] != EMPTY) h= (h + 1) & MASK;
   u32 fp= u32(cur >> BITS) & FP_MASK;
   tab[h]= (fp << VBITS) | (k + 1);
   cur= mul(cur, G_65537);
  }
 }
 u32 lookup(u64 t) const {
  u32 h= u32(t) & MASK;
  u32 want= (u32(t >> BITS) & FP_MASK) << VBITS;
  while(tab[h] != EMPTY) {
   u32 e= tab[h];
   if(((e ^ want) >> VBITS) == 0) return (e & VMASK) - 1;
   h= (h + 1) & MASK;
  }
  return ~u32(0);
 }
};
inline DirectLogTable_65537 direct_65537;
// =============================================================================
// BSGSTable6700417: 6700417 元 subgroup 用 BSGS (packed layout).
//
// レイアウト: u64 entry = (key & HI_MASK) | value
//   - 下位 17 bit: value (j ∈ [0, m))
//   - bits 17..18: gap (常に 0)
//   - 上位 45 bit: key の fingerprint (= key & ~MASK)
//
// 比較は ((e ^ t) & HI_MASK) == 0 を、コンパイラは自動的に (e ^ t) <= MASK
// に書き換える (5 命令のホットループになる)。
//
// inv_base_m = G_6700417^(P_BIG - m) は事前計算した constexpr 定数を使用。
// =============================================================================
struct BSGSTable6700417 {
 static constexpr u32 mask= 524287;  // 19-bit, cap = 524288
 static constexpr u64 q= P_BIG;
 static constexpr u32 m= 131072;
 static constexpr u64 max_i= 52;
 static constexpr u64 HI_MASK= ~u64(mask);  // = 0xFFFFFFFFFFF80000
 static constexpr u64 EMPTY= ~u64(0);
 static constexpr u64 inv_base_m= 0x1489880b9cf723deull;
 std::vector<u64> tab;
 void build() {
  u64 base= G_6700417;
  tab.assign(mask + 1, EMPTY);
  u64 cur= 1;
  for(u32 j= 0; j < m; ++j) {
   u32 h= u32(cur) & mask;
   while(tab[h] != EMPTY) h= (h + 1) & mask;
   tab[h]= (cur & HI_MASK) | u64(j);
   cur= mul(cur, base);
  }
 }
 u32 solve(u64 target) const {
  u64 t= target;
  for(u8 i= 0; i <= max_i; ++i) {
   u32 h= u32(t) & mask;
   while(tab[h] != EMPTY) {
    u64 e= tab[h];
    if(((e ^ t) & HI_MASK) == 0) return u32(i * m + u32(e & mask));
    h= (h + 1) & mask;
   }
   t= mul(t, inv_base_m);
  }
  return q;  // not found
 }
};
inline BSGSTable6700417 bsgs_6700417;

inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 // PerfectHash641 は constexpr で .rodata に焼かれているため build 不要。
 direct_65537.build();
 bsgs_6700417.build();
}
u64 log_g(u64 x) {
 assert(x);
 u64 N, s, x_f16, x_65537, x_6700417, x_641;
 __m256i Ns= mul2(_mm256_set_epi64x(0, sq(x), 0, frob32(x)), _mm256_set1_epi64x(x), N, s);
 mul2(_mm256_set_epi64x(0, frob2(s), 0, frob16(N)), Ns, x_f16, s);
 s= mul(s, frob4(s));
 s= mul(s, frob8(s));  // 2^16-1
 mul2(_mm256_set_epi64x(0, frob16(s), 0, frob32(s)), _mm256_set1_epi64x(s), x_65537, s);
 u64 T2= sq(s), T3= mul(s, T2), T5, T10;
 u64 acc= T3;
 mul2(_mm256_set_epi64x(0, T2, 0, frob4(acc)), _mm256_set1_epi64x(T3), acc, T5);
 mul2(_mm256_set_epi64x(0, T5, 0, frob7(acc)), _mm256_set1_epi64x(T5), acc, T10);
 acc= mul(acc, T10);
 acc= mul(frob3(acc), T3);
 mul2(_mm256_set_epi64x(0, frob7(acc), 0, frob7(T5)), _mm256_set1_epi64x(s), x_6700417, x_641);
 const u16 r1= LN16[u16(x_f16)];
 const u32 r0= PerfectHash641::lookup(x_641);
 const u32 r2= direct_65537.lookup(x_65537);
 const u32 r3= bsgs_6700417.solve(x_6700417);
 const u16 cur_mod_641= r1 % P_641;
 const u16 diff0= (r0 + P_641 - cur_mod_641) % P_641;
 const u16 t0= (diff0 * INV_65535_641) % P_641;
 const u32 cur_after0_mod_F17= (r1 + MOD_F16 * t0) % P_F17;
 const u32 diff2= (r2 + P_F17 - cur_after0_mod_F17) % P_F17;
 const u32 t2= (diff2 * INV_MOD2_F17) % P_F17;
 const u32 cur_after2_mod_BIG= (r1 + (MOD_F16 * t0) % P_BIG + (MOD_F16_641 * t2) % P_BIG) % P_BIG;
 const u32 diff3= (r3 + P_BIG - cur_after2_mod_BIG) % P_BIG;
 const u32 t3= (u64(diff3) * INV_MOD3_BIG) % P_BIG;
 return u64(r1) + MOD_F16 * t0 + MOD_F16_641 * t2 + MOD_F16_641_F17 * t3;
}
}  // namespace gf2_64_log_pohlig_v11
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& xs) {
  using gf2_64_log_pohlig_v11::init_tables;
  using gf2_64_log_pohlig_v11::log_g;
  init_tables();
  vector<u64> ans(xs.size());
  for(size_t i= 0; i < xs.size(); ++i) ans[i]= log_g(xs[i]);
  return ans;
 }
};
