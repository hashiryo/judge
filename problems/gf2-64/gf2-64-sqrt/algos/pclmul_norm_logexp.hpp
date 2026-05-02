#pragma once
// Nimber 流の sqrt を poly basis で再現したバリアント。
// (linear_map.hpp が圧倒的に速いはずだが、比較用として実装)
//
// Nimber の sqrt:
//   入力 α = a_0 + a_1 ω + a_2 ω^2 + a_3 ω^3 (a_i ∈ F_{2^16}, nim 基底)
//   1) half 補正: a_1 ^= half(a_3 ^ a_2),  a_2 ^= half(a_3),
//                 a_0 ^= half(a_1) ^ half<6>(a_3)
//   2) 各成分に F_{2^16} 平方根: sqrt_f16(a_i) = pw[(65537 * ln[a_i]) >> 1]
//   3) 再合成
//
// 我々の F_{2^16} log/exp テーブルを使って同じ計算を行う。basis change は
// pclmul_norm_logexp.hpp と同じ仕組み。
#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("pclmul,bmi2")
#endif
#include "_common.hpp"
#include "../../_shared/pclmul_core.hpp"
#include "../../_shared/basis_change.hpp"

#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#define PCLMUL_RUN [[gnu::target("pclmul,bmi2")]]
#else
#define PCLMUL_RUN
#endif

namespace gf2_64_sqrt_norm_logexp {
using u16 = unsigned short;

inline u16 PW16[65536], LN16[65536];
inline bool inited = false;

PCLMUL_FN void init_tables() {
 if (inited) return;
 inited = true;
 // F_{2^16} log/exp (Nimber.hpp 互換)
 PW16[0] = PW16[65535] = 1;
 for (int i = 1; i < 65535; ++i) {
  PW16[i] = u16((PW16[i-1] << 1) ^ (0x1681fu & u16(-(PW16[i-1] >= 0x8000u))));
 }
 constexpr u16 f2n[16] = {0x0001u, 0x2827u, 0x392bu, 0x8000u, 0x20fdu, 0x4d1du, 0xde4au, 0x0a17u,
                          0x3464u, 0xe3a9u, 0x6d8du, 0x34bcu, 0xa921u, 0xa173u, 0x0ebcu, 0x0e69u};
 for (int i = 1; i < 65535; ++i) {
  u16 x = PW16[i], y = 0;
  for (; x; x &= x - 1) y ^= f2n[__builtin_ctz(x)];
  PW16[i] = y;
  LN16[y] = u16(i);
 }
 LN16[1] = 0;
}

// half<h>(A) = pw[(ln[A] + h) % 65535]   (Nimber では half の h は塔の補正定数)
template<u16 h>
[[gnu::always_inline]] inline u16 half(u16 A) {
 return A ? PW16[(u32(LN16[A]) + h) % 65535] : 0;
}

// F_{2^16} sqrt: A^{2^15} を log/exp で
[[gnu::always_inline]] inline u16 sqrt_f16(u16 A) {
 return A ? PW16[u16((65537u * u32(LN16[A])) >> 1)] : 0;
}

PCLMUL_FN u64 sqrt_via_tower(u64 a_poly) {
 const u64 a_nim = gf2_64_basis::poly_to_nim(a_poly);
 u16 a0 = u16(a_nim), a1 = u16(a_nim >> 16), a2 = u16(a_nim >> 32), a3 = u16(a_nim >> 48);
 // half 補正 (Nimber.hpp と同じ順序・定数)
 a1 ^= half<3>(u16(a3 ^ a2));
 a2 ^= half<3>(a3);
 a0 ^= half<3>(a1) ^ half<6>(a3);
 // 各成分に F_{2^16} sqrt
 const u16 b0 = sqrt_f16(a0), b1 = sqrt_f16(a1), b2 = sqrt_f16(a2), b3 = sqrt_f16(a3);
 const u64 b_nim = u64(b0) | (u64(b1) << 16) | (u64(b2) << 32) | (u64(b3) << 48);
 return gf2_64_basis::nim_to_poly(b_nim);
}

} // namespace gf2_64_sqrt_norm_logexp

struct GF2_64Op {
 PCLMUL_RUN static vector<u64> run(const vector<u64>& as) {
  using gf2_64_sqrt_norm_logexp::sqrt_via_tower;
  using gf2_64_sqrt_norm_logexp::init_tables;
  init_tables();
  vector<u64> ans(as.size());
  for (size_t i = 0; i < as.size(); ++i) ans[i] = sqrt_via_tower(as[i]);
  return ans;
 }
};
