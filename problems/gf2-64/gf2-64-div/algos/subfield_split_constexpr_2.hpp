// INV_LOW[low16(a)] = low16(a^{-1}) を constexpr で構築
// 元のコードの「u16(mul(g, a)) で直接引く」スタイルと互換

#include "../../_shared/_common.hpp"
#include "../../_shared/sq.hpp"
#include "../../_shared/frob.hpp"

using gf2_64_pclmul::frob16;
using gf2_64_pclmul::frob32;
using gf2_64_pclmul::frob4;
using gf2_64_pclmul::frob48;
using gf2_64_pclmul::mul;
using gf2_64_pclmul::sq;
struct LinTable {
 u16 T_lo[256];
 u16 T_hi[256];
};
// 行列 M (16x16) を 2段ルックアップに変換
// (M を「列」表現として持ち、入力ビットが1の列を XOR)
constexpr LinTable mat_to_table(const array<u16, 16>& M) {
 // M の i 列 = (M.row[r].bit(i) for r=0..15) を 16bit にまとめる
 u16 col[16]= {};
 for(int i= 0; i < 16; ++i)
  for(int r= 0; r < 16; ++r)
   if((M[r] >> i) & 1) col[i]|= 1 << r;

 LinTable t{};
 for(int v= 0; v < 256; ++v) {
  u16 lo= 0, hi= 0;
  for(int j= 0; j < 8; ++j) {
   if((v >> j) & 1) {
    lo^= col[j];
    hi^= col[j + 8];
   }
  }
  t.T_lo[v]= lo;
  t.T_hi[v]= hi;
 }
 return t;
}
// === メイン: low16 → low16 の inverse テーブル ===
constexpr auto INV_LOW= []() {
 constexpr array<u16, 16> M_NAT_TO_LOW= {59577U, 18906U, 52628U, 23944U, 45524U, 3436U, 57636U, 59944U, 41644U, 45342U, 50170U, 27798U, 34196U, 59186U, 37224U, 53912U};
 constexpr auto NAT_TO_LOW= mat_to_table(M_NAT_TO_LOW);
 array<u16, 65536> t{};
 // σ chain を natural 表現で 1 周だけ走らせ nat[k] = σ^k を保持。
 // log/exp テーブル経由を廃止 → INV_LOW を nat 空間でじか埋めする。
 u16 nat[65535];
 {
  u16 cur= 1;
  for(uint32_t i= 0; i < 65535; ++i) {
   nat[i]= cur;
   cur= (u16)(u16(cur << 1) ^ u16(0x002DU & -u16(cur >> 15)));
  }
 }
 // σ^0 = 1: INV_LOW[low(1)] = low(1)
 u16 lo1= NAT_TO_LOW.T_lo[1] ^ NAT_TO_LOW.T_hi[0];
 t[lo1]= lo1;
 // pair (k, 65535-k) で 2 entry ずつ書く → 反復は 32767 回で済む
 for(uint32_t k= 1; k <= 32767; ++k) {
  u16 nk= nat[k];
  u16 nik= nat[65535 - k];
  u16 lo_k= NAT_TO_LOW.T_lo[u8(nk)] ^ NAT_TO_LOW.T_hi[nk >> 8];
  u16 lo_ik= NAT_TO_LOW.T_lo[u8(nik)] ^ NAT_TO_LOW.T_hi[nik >> 8];
  t[lo_k]= lo_ik;
  t[lo_ik]= lo_k;
 }
 return t;
}();
// === embed テーブル(low16 → 64bit、SUBFIELD_BASIS による)===
constexpr array<u64, 16> SUBFIELD_BASIS= {1ULL, 6899425322512154626ULL, 12712641506861907972ULL, 12687683756412895240ULL, 13108774640850436112ULL, 1196746230653255712ULL, 13779846473293824064ULL, 1136705091741089920ULL, 13132935623751303424ULL, 12256911237861802496ULL, 1968662052679910400ULL, 13476734309037115392ULL, 31478309824172032ULL, 5397840376063860736ULL, 18145356609018085376ULL, 2133828226494464000ULL};
constexpr auto EMBED= []() {
 array<array<u64, 256>, 2> t{};
 for(int half= 0; half < 2; ++half) {
  for(int i= 0; i < 256; ++i) {
   u64 v= 0;
   for(int b= 0; b < 8; ++b)
    if((i >> b) & 1) v^= SUBFIELD_BASIS[b + half * 8];
   t[half][i]= v;
  }
 }
 return t;
}();
constexpr u64 embed_idx(u16 idx) { return EMBED[0][u8(idx)] ^ EMBED[1][idx >> 8]; }
u64 inv(u64 a) {
 assert(a != 0);
 u64 g= mul(frob16(a), mul(frob32(a), frob48(a)));
 u64 b= embed_idx(INV_LOW[u16(mul(g, a))]);
 return mul(b, g);
}
struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_pclmul::mul;
  vector<u64> ans(as.size());
  for(size_t i= 0; i < as.size(); ++i) ans[i]= mul(as[i], inv(bs[i]));
  return ans;
 }
};
