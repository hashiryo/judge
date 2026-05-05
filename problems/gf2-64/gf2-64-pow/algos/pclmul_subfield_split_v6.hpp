#pragma once
// v5_5 の発展形: G_b (order 65537) の hash log を撤廃し、α^e の分解を二重に行う。
//
// 第 1 段 (v4_4 / v5_5 と同じ):
//   α^e = N2(α)^q · α^r   (q = e / M_VAL2, r = e mod M_VAL2, M_VAL2 = 2^32 + 1)
//   N2(α) = α · α^{2^32} ∈ F_{2^32}^* (order 65535·65537), q ∈ [0, 2^32 - 1]
// 第 2 段 (v6 で追加):
//   q = q_2 · (2^16 + 1) + r_2  (q_2 ∈ [0, 65535], r_2 ∈ [0, 65537))
//   N2(α)^q = N2(α)^{(2^16+1)·q_2} · N2(α)^{r_2}
//           = N(α)^{q_2}              · N2(α)^{r_2}
//   N(α) = N2(α)^{2^16+1} = α · α^{2^16} · α^{2^32} · α^{2^48} ∈ F_{2^16}^*
//
//   - N(α)^{q_2} : 既存 LN_SIGMA / PW_SIGMA_IDX で 1 lookup + 1 mul + 1 lookup (~6 cycle)
//   - N2(α)^{r_2}: r_2 が高々 17 bit なので pow_byte_window で計算 (~80 cycle)
//
// 利点: G_b 用の hash log table (1 MiB) と PW_H2 (512 KiB) を完全撤廃 → メモリ 1.5 MiB 削減。
// 欠点: per-query で N2^{r_2} の pow_byte_window を 1 回追加 (~80 cycle)。
#pragma GCC optimize("O3,unroll-loops")
#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"
namespace gf2_64_pow_subfield_split_v6 {
using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;

// constexpr で [2][256] EMBED_BYTE テーブルを構築 (v5_5 と同一)
constexpr auto EMBED_BYTE= []() {
 constexpr std::array<uint64_t, 16> SUBFIELD_BASIS= {0x0000000000000001ULL, 0x5fbfaec6aeac0002ULL, 0xb06c601895640004ULL, 0xb013b5277b7c0008ULL, 0xb5ebb915248a0010ULL, 0x109bb25b2c600020ULL, 0xbf3bd95bd4190040ULL, 0x0fc66342279b0080ULL, 0xb6418f5e57c50100ULL, 0xaa194bd4b83f0200ULL, 0x1b5217b4dcc70400ULL, 0xbb06fa73867a0800ULL, 0x006fd55b23331000ULL, 0x4ae8fb39198c2000ULL, 0xfbd141b29b4f4000ULL, 0x1d9ce1776be78000ULL};
 std::array<std::array<uint64_t, 256>, 2> T{};
 for(int half= 0; half < 2; ++half) {
  for(int i= 0; i < 256; ++i) {
   uint64_t v= 0;
   for(int b= 0; b < 8; ++b) {
    if((i >> b) & 1) v^= SUBFIELD_BASIS[b + half * 8];
   }
   T[half][i]= v;
  }
 }
 return T;
}();
inline u64 embed_idx(u16 idx) { return EMBED_BYTE[0][u8(idx)] ^ EMBED_BYTE[1][u8(idx >> 8)]; }
constexpr u64 M_VAL2= (u64(1) << 32) + 1;
constexpr u32 M_F17= 65537;  // 2^16 + 1

inline u16 LN_SIGMA[65536];      // LN_SIGMA[u16(σ^k)] = k  (raw log)
inline u16 PW_SIGMA_IDX[65536];  // PW_SIGMA_IDX[k] = u16(σ^k)
u64 pow_byte_window(u64 g, u64 e) {
 if(!e) return 1;
 u64 T[16]= {1, g};
 for(int i= 2; i < 16; ++i) T[i]= mul(T[i - 1], g);
 int top= 15 - (__builtin_clzll(e) >> 2);
 u64 acc= T[(e >> (4 * top)) & 0xF];
 for(int i= top - 1; i >= 0; --i) {
  acc= frob4(acc);
  u8 chunk= (e >> (4 * i)) & 0xF;
  if(chunk) acc= mul(acc, T[chunk]);
 }
 return acc;
}
inline bool inited= false;
void init_tables() {
 if(inited) return;
 inited= true;
 constexpr u64 SIGMA= 0xa1573a4da2bc3a32ull;
 u64 cur= 1;
 for(u16 k= 0; k < 65535; ++k) {
  u16 idx= u16(cur);
  LN_SIGMA[idx]= k;
  PW_SIGMA_IDX[k]= idx;
  cur= mul(cur, SIGMA);
 }
 LN_SIGMA[0]= 0;
 PW_SIGMA_IDX[65535]= 1;
}
u64 pow(u64 a, u64 e) {
 if(!e) return 1;
 if(!a) return 0;
 const u32 q= u32(e / M_VAL2);
 if(!q) return pow_byte_window(a, e);
 const u64 r= e - u64(q) * M_VAL2;
 // q = q_2 · (2^16+1) + r_2
 const u32 q_2= q / M_F17;           // ∈ [0, 65535]
 const u32 r_2= q - q_2 * M_F17;     // ∈ [0, 65537)
 const u64 N2= mul(a, frob32(a));    // N2(α) ∈ F_{2^32}^*
 const u64 Na= mul(N2, frob16(N2));  // N(α) = N2(α)^{2^16+1} ∈ F_{2^16}^*
 // b1 = N(α)^{q_2}  (q_2 = 0 でも LN_SIGMA[Na]·0 = 0 → PW_SIGMA_IDX[0] = 1 で正しく)
 const u32 L_a= LN_SIGMA[u16(Na)];
 const u32 e_a= u64(q_2) * L_a % 65535u;
 const u64 b1= embed_idx(PW_SIGMA_IDX[e_a]);
 // b2 = N2(α)^{r_2}  (r_2 ≤ 65537、~17 bit)
 const u64 b2= pow_byte_window(N2, r_2);
 const u64 b= mul(b1, b2);
 if(!r) return b;
 const u64 g= pow_byte_window(a, r);
 return mul(b, g);
}
}  // namespace gf2_64_pow_subfield_split_v6
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& es) {
  using gf2_64_pow_subfield_split_v6::init_tables;
  using gf2_64_pow_subfield_split_v6::pow;
  init_tables();
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= pow(as[i], es[i]);
  return ans;
 }
};
