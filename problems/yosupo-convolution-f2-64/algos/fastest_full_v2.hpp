#pragma once
#include "_common.hpp"
// F_{2^64} 上の additive (nim) FFT に対する全部入り最適化版。
// アルゴリズム: Lin-Chung-Han / Cantor 系の characteristic-2 FFT (Phase A: 線形
// XOR-only butterfly + 並べ替え、 Phase B: twiddle 乗算)。
//
// 取り入れた最適化 (積み重ね):
//
// 1. Cantor 対称性 (基底の末尾 = 1 → twiddle 表の上下半分が +1 差)
//      a[i][j + 2^(i-1)] = a[i][j] ^ 1
//    これにより butterfly 1 個あたり 2 mul → 1 mul に半減:
//      hi * g0[k]           (= p)
//      hi * g1[k] = p ^ hi  (∵ GF(2) 上の分配律)
//    IFFT 側も同様: f0_new = (lo+hi)*g0 + lo, f1_new = lo+hi で 1 mul に削減。
//
// 2. halftable: 上半分は g0[k] ^ 1 で復元できるので保存不要。 twiddle メモリ半減
//    (2^(n+1) - 1 → 2^n - 1 entry, s=2^20 で 16MB → 8MB)、 L3 圏内に収まりやすく
//    memory-bandwidth 改善。
//
// 3. g0[0] = 0 スキップ + mul2 ペアリング: 各ブロック先頭の i=0 は PCLMUL 不要、
//    残り (i=1..half-1) を 2 個ずつ VPCLMULQDQ で並列化。
//
// 4. butterfly pair の SIMD pipeline (btf_fft_pair / btf_ifft_pair):
//    load → 256-bit 展開 → VPCLMULQDQ + 並列 reduction → __m256i のまま XOR →
//    pack して __m128i 一括 store。 vector ↔ scalar の往復を排除し、 ペアあたり
//    extract / scalar XOR / scalar store を ~7-8 命令削減。
//
// 5. Phase A の最終 2 段 + shuffle 融合: 連続する 3 パス (m=2 stage, m=1 stage,
//    shuffle) は全部 F_2 線形変換なので、 8-tuple に対する 1 パスに圧縮。
//    メモリトラフィックが ~3x 減 (36 ops → 16 ops / 8-tuple)。 IFFT 側も対称に。
//    len=2 の Phase A iteration は identity + swap で実質 no-op なので skip。
//
// 6. init() の per-level 動的計算を排除: descent 中の basis が常に CHAIN suffix
//    と一致 (`b.back() = c[62] = 1` が不変)、 inv/mul/sq/pop は全部 no-op だった。
//    constexpr で焼いた CHAIN を直接参照する Gray code build だけに簡略化。
//
// 必要拡張: PCLMUL + VPCLMULQDQ + AVX2 + BMI2 (Ice Lake / Zen3 以降)。

#pragma GCC optimize("O3,unroll-loops")
#include "../../gf2-64/_shared/_common.hpp"
#include "../../gf2-64/_shared/mul.hpp"
namespace conv_f2_64_full_v2 {

using gf2_64_pclmul::mul;
// =============================================================================
// constexpr GF(2^64) mul / sq (chain を .rodata に焼くため)
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
// Artin-Schreier 連鎖 c[k] = P^k(2), P(x) = x²+x。
// 実測で c[62] = 1, c[63] = 0 なので非零 basis は c[0..62] の 63 個。
// =============================================================================
constexpr int CHAIN_LEN= 63;
constexpr auto CHAIN= []() {
 array<u64, CHAIN_LEN> c{};
 c[0]= 2;
 for(int k= 1; k < CHAIN_LEN; ++k) c[k]= sq_ce(c[k - 1]) ^ c[k - 1];
 return c;
}();
// =============================================================================
// fused butterfly pair (FFT / IFFT) — load → expand → mul + reduction → XOR → store
// を __m256i のまま実行。 vector ↔ scalar 往復を排除。
// =============================================================================
//
// __m128i (a0, a1) → __m256i (a0, _, a1, _) への展開 (vpermq 1 命令)
// _MM_SHUFFLE(d, c, b, a) は output[k] = src[k番目引数] を意味し、 ここでは
// output = (src[0], src[2]=0, src[1], src[3]=0) が欲しいので imm = _MM_SHUFFLE(3, 1, 2, 0)。
VPCLMUL inline __m256i expand_pair(const u64* p) {
 __m128i v= _mm_loadu_si128((const __m128i*)p);
 return _mm256_permute4x64_epi64(_mm256_castsi128_si256(v), _MM_SHUFFLE(3, 1, 2, 0));
}
// __m256i (x0, _, x1, _) → __m128i (x0, x1) (lane 0, 2 を low 128 へ集約)
// imm = _MM_SHUFFLE(3, 2, 2, 0): output = (src[0], src[2], src[2], src[3])。 low128 = (src[0], src[2])。
VPCLMUL inline __m128i pack_pair(__m256i v) { return _mm256_castsi256_si128(_mm256_permute4x64_epi64(v, _MM_SHUFFLE(3, 2, 2, 0))); }
// VPCLMULQDQ + 並列 reduction を __m256i で完結 (mul2 と同 idiom、 結果は (p0, _, p1, _))

inline const __m256i RED_TABLE= _mm256_setr_epi8(0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 27, 45, 54, 90, 65, 119, 108, 0, 0, 0, 0, 0, 0, 0, 0);
VPCLMUL inline __m256i clmul_reduce_pair(__m256i a_vec, __m256i b_vec) {
 __m256i prod= _mm256_clmulepi64_epi128(a_vec, b_vec, 0);
 __m256i d_full= _mm256_xor_si256(prod, _mm256_slli_epi64(prod, 1));
 __m256i red1_full= _mm256_xor_si256(d_full, _mm256_slli_epi64(d_full, 3));
 __m256i red1_shift= _mm256_srli_si256(red1_full, 8);
 __m256i h_idx= _mm256_srli_epi64(prod, 60);
 __m256i indices= _mm256_srli_si256(h_idx, 8);
 __m256i red_vec= _mm256_shuffle_epi8(RED_TABLE, indices);
 return _mm256_xor_si256(_mm256_xor_si256(prod, red1_shift), red_vec);
}
// FFT 用 butterfly pair: (i, i+1) を一括処理。
// 計算: t = lo ^ hi*g0; f0' = t; f1' = t ^ hi
VPCLMUL inline void btf_fft_pair(u64* f0_ptr, u64* f1_ptr, const u64* g0_ptr) {
 __m256i hi_vec= expand_pair(f1_ptr);
 __m256i g_vec= expand_pair(g0_ptr);
 __m256i p_vec= clmul_reduce_pair(hi_vec, g_vec);
 __m256i lo_vec= expand_pair(f0_ptr);
 __m256i t_vec= _mm256_xor_si256(lo_vec, p_vec);
 __m256i f1_new= _mm256_xor_si256(t_vec, hi_vec);
 _mm_storeu_si128((__m128i*)f0_ptr, pack_pair(t_vec));
 _mm_storeu_si128((__m128i*)f1_ptr, pack_pair(f1_new));
}
// IFFT 用 butterfly pair:
// 計算: s = lo ^ hi; f0' = lo ^ s*g0; f1' = s
VPCLMUL inline void btf_ifft_pair(u64* f0_ptr, u64* f1_ptr, const u64* g0_ptr) {
 __m256i lo_vec= expand_pair(f0_ptr);
 __m256i hi_vec= expand_pair(f1_ptr);
 __m256i s_vec= _mm256_xor_si256(lo_vec, hi_vec);
 __m256i g_vec= expand_pair(g0_ptr);
 __m256i p_vec= clmul_reduce_pair(s_vec, g_vec);
 __m256i f0_new= _mm256_xor_si256(lo_vec, p_vec);
 _mm_storeu_si128((__m128i*)f0_ptr, pack_pair(f0_new));
 _mm_storeu_si128((__m128i*)f1_ptr, pack_pair(s_vec));
}
template <class T> int msb(T n) { return n == 0 ? -1 : 63 - __builtin_clzll(n); }
template <class T> T ceil_pow2(T n) { return n <= 1 ? T(1) : T(1) << (msb(n - 1) + 1); }
// =============================================================================
// Twiddle: 下半分のみ (上半分は g[i]^1 = g[i+half] で復元可能なので不要)
//   level i (i ≥ 1) の table 長 = 2^(i-1)
//   level 0 はサイズ 1 の唯一の要素 {0} (ホットパスから読まれない)
// =============================================================================
struct nim_fft_data {
 std::vector<std::vector<u64>> a;
 void init(int n) {
  if((int)a.size() > n) return;
  a.resize(n + 1);
  // level i での basis = CHAIN[CHAIN_LEN-i .. CHAIN_LEN-1]。
  // back = CHAIN[CHAIN_LEN-1] = c[62] = 1 (Cantor 正規化済み) なので、
  // 旧コードの per-level inv / mul / sq / pop は全部 no-op。 直接参照で十分。
  for(int i= 0; i <= n; ++i) {
   auto& na= a[i];
   const int sz= (i == 0) ? 1 : (1 << (i - 1));
   na.resize(sz);
   na[0]= 0;
   const u64* basis= &CHAIN[CHAIN_LEN - i];  // basis[k] = level-i basis の k 番目
   for(int j= 1; j < sz; ++j) {
    int k= __builtin_ctz(j);
    na[j]= na[j ^ (1 << k)] ^ basis[k];
   }
  }
 }
};
inline nim_fft_data nim_data;
// =============================================================================
// nim FFT (u64 そのまま、 + は ^、 * は mul / mul2)
// =============================================================================
VPCLMUL inline void nim_fft(std::vector<u64>& f) {
 int n= (int)f.size();
 std::vector<u64> f2(n);
 int len= n;
 // Phase A: 最終 2 段 + shuffle 融合。 早期段 (m ≥ 4) は inplace、 最後 m=2,m=1
 // と shuffle は 8-tuple 単位で 1 パスに集約。
 while(len > 1) {
  if(len == 2) {
   // 旧コードでは「m=0 の butterfly skip + identity shuffle + swap」 で f に
   // 元の内容が残るだけの no-op (1 pass のメモリ無駄)。 まるごとスキップ。
   len/= 2;
   continue;
  }
  if(len == 4) {
   // 1 stage (m=1) + shuffle 融合: 4-tuple p[0..3] → out
   // m=1: (a, b^c^d, c^d, d). shuffle: out = [a, c^d, b^c^d, d]
   for(int l= 0; l < n; l+= 4) {
    u64 p0= f[l], p1= f[l + 1], p2= f[l + 2], p3= f[l + 3];
    u64 t= p2 ^ p3;
    f2[l]= p0;
    f2[l + 1]= t;
    f2[l + 2]= p1 ^ t;
    f2[l + 3]= p3;
   }
   std::swap(f, f2);
   len/= 2;
   continue;
  }
  // len >= 8: 早期段 m = len/4, ..., 4 を inplace。 m=2 と m=1 と shuffle は融合。
  for(int l= 0; l < n; l+= len) {
   for(int m= len / 4; m >= 4; m/= 2) {
    for(int t= 0; t < len; t+= m * 4) {
     for(int i= 0; i < m; ++i) {
      u64 b= f[l + t + m + i], c= f[l + t + m * 2 + i], d= f[l + t + m * 3 + i];
      f[l + t + m + i]= b ^ c ^ d;
      f[l + t + m * 2 + i]= c ^ d;
     }
    }
   }
   const int hf= len / 2;
   // 8-tuple s = 0..len/8-1, src f[l + 8s..8s+7], dst f2[l + 4s + {0..3, hf..hf+3}]
   for(int s= 0; s < len / 8; ++s) {
    int b8= l + s * 8;
    u64 p0= f[b8], p1= f[b8 + 1], p2= f[b8 + 2], p3= f[b8 + 3];
    u64 p4= f[b8 + 4], p5= f[b8 + 5], p6= f[b8 + 6], p7= f[b8 + 7];
    u64 a3457= p3 ^ p4 ^ p5;         // 共有: out[1] と out[4] のかたまり
    u64 A2_7= p2 ^ a3457 ^ p6 ^ p7;  // p2^p3^p4^p5^p6^p7
    int d0= l + s * 4;
    int d1= d0 + hf;
    f2[d0 + 0]= p0;
    f2[d0 + 1]= A2_7;
    f2[d0 + 2]= p4 ^ p6;
    f2[d0 + 3]= p6 ^ p7;
    f2[d1 + 0]= p1 ^ A2_7;
    f2[d1 + 1]= p3 ^ p5 ^ p7;
    f2[d1 + 2]= p5 ^ p6;
    f2[d1 + 3]= p7;
   }
  }
  std::swap(f, f2);
  len/= 2;
 }
 // Phase B: Cantor 対称性 (g[i+half]=g[i]^1) + g0[0]=0 スキップ + mul2 ペアリング
 while(len < n) {
  len*= 2;
  const std::vector<u64>& g= nim_data.a[msb(len)];
  for(int l= 0; l < n; l+= len) {
   const int half= len / 2;
   const u64* g0= g.data();
   u64* f0= f.data() + l;
   u64* f1= f.data() + l + half;
   // i = 0: g0[0] = 0 → PCLMUL も g0 ロードも不要
   {
    u64 lo= f0[0], hi= f1[0];
    f0[0]= lo;
    f1[0]= lo ^ hi;
   }
   // i = 1..half-1 を 2 個ずつ pair で fused SIMD butterfly
   int i= 1;
   for(; i + 1 < half; i+= 2) {
    btf_fft_pair(f0 + i, f1 + i, g0 + i);
   }
   // tail: 残った 1 個 (half が偶数なら必ず i = half-1 が単独で残る)
   for(; i < half; ++i) {
    u64 lo= f0[i], hi= f1[i];
    u64 p= mul(hi, g0[i]);
    u64 t= lo ^ p;
    f0[i]= t;
    f1[i]= t ^ hi;
   }
  }
 }
}
VPCLMUL inline void nim_ifft(std::vector<u64>& f) {
 int n= (int)f.size();
 std::vector<u64> f2(n);
 int len= n;
 // IFFT: Cantor 対称性 + g0[0]=0 スキップ + mul2 ペアリング
 while(len > 1) {
  const std::vector<u64>& g= nim_data.a[msb(len)];
  for(int l= 0; l < n; l+= len) {
   const int half= len / 2;
   const u64* g0= g.data();
   u64* f0= f.data() + l;
   u64* f1= f.data() + l + half;
   // i = 0: g0[0] = 0 → s*g0 = 0, f0 = lo, f1 = s = lo^hi
   {
    u64 lo= f0[0], hi= f1[0];
    f0[0]= lo;
    f1[0]= lo ^ hi;
   }
   int i= 1;
   for(; i + 1 < half; i+= 2) {
    btf_ifft_pair(f0 + i, f1 + i, g0 + i);
   }
   for(; i < half; ++i) {
    u64 lo= f0[i], hi= f1[i];
    u64 s= lo ^ hi;
    u64 p= mul(s, g0[i]);
    f0[i]= lo ^ p;
    f1[i]= s;
   }
  }
  len/= 2;
 }
 // IFFT Phase A: shuffle + 最初 2 段 (m=1, m=2) を 8-tuple 単位で融合。
 // m ≥ 4 の段は融合後に inplace で実行。
 while(len < n) {
  len*= 2;
  if(len == 2) {
   // shuffle (identity) + butterfly (m=0、 入らない) なので no-op
   continue;
  }
  if(len == 4) {
   // shuffle + m=1 融合: q[0]=src[0], q[1]=src[2], q[2]=src[1], q[3]=src[3]
   //   out[0]=q[0], out[1]=q[1]^q[2], out[2]=q[2]^q[3], out[3]=q[3]
   for(int l= 0; l < n; l+= 4) {
    u64 s0= f[l], s1= f[l + 1], s2= f[l + 2], s3= f[l + 3];
    f2[l]= s0;
    f2[l + 1]= s2 ^ s1;
    f2[l + 2]= s1 ^ s3;
    f2[l + 3]= s3;
   }
   std::swap(f, f2);
   continue;
  }
  // len >= 8: shuffle + m=1 + m=2 を 8-tuple 単位で融合 → f2、 swap、 m ≥ 4 を inplace
  const int hf= len / 2;
  for(int l= 0; l < n; l+= len) {
   for(int s= 0; s < len / 8; ++s) {
    int b4lo= l + s * 4;
    int b4hi= b4lo + hf;
    u64 q0= f[b4lo + 0], q2= f[b4lo + 1], q4= f[b4lo + 2], q6= f[b4lo + 3];
    u64 q1= f[b4hi + 0], q3= f[b4hi + 1], q5= f[b4hi + 2], q7= f[b4hi + 3];
    u64 c67= q6 ^ q7;
    u64 c56= q5 ^ q6;
    int d8= l + s * 8;
    f2[d8 + 0]= q0;
    f2[d8 + 1]= q1 ^ q2;
    f2[d8 + 2]= q2 ^ q3 ^ q4;
    f2[d8 + 3]= q3 ^ c56;
    f2[d8 + 4]= q4 ^ c67;
    f2[d8 + 5]= q5 ^ c67;
    f2[d8 + 6]= c67;
    f2[d8 + 7]= q7;
   }
  }
  std::swap(f, f2);
  // 残りの段 m = 4, 8, ..., len/4 を inplace
  for(int l= 0; l < n; l+= len) {
   for(int m= 4; m <= len / 4; m*= 2) {
    for(int t= 0; t < len; t+= m * 4) {
     for(int i= 0; i < m; ++i) {
      u64 b= f[l + t + m + i], c= f[l + t + m * 2 + i], d= f[l + t + m * 3 + i];
      f[l + t + m + i]= b ^ c;
      f[l + t + m * 2 + i]= c ^ d;
     }
    }
   }
  }
 }
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
}  // namespace conv_f2_64_full_v2
struct Solver {
 VPCLMUL static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_full_v2;
  auto c= nim_convolution(a_in, b_in);
  return c;
 }
};
