// Special judge for gf2-64-log.
// 入力: T 個の x_i
// 期待出力 (コンテスト想定): 各 k_i with g^k_i = x_i (g = 2 固定、原始元)
// ジャッジ: ユーザの k_user に対し g^{k_user} == x_i であれば AC。
//
// (g は P(x) = x^64+x^4+x^3+x+1 の下で原始元なので k は [0, 2^64-2] で一意。
//  実装毎に同値表現が違う可能性に備えて、checker で再計算検証する。)
#include "testlib.h"
#include <vector>

using u64 = unsigned long long;
constexpr u64 IRRED_LOW = 0x1Bu;
constexpr u64 G = 2;

// 素朴 GF(2^64) 演算 (testlib 環境で __m128i 使えない可能性に備えて pure scalar)。
static inline std::pair<u64, u64> clmul_loop(u64 a, u64 b) {
 u64 lo = 0, hi = 0;
 for (int i = 0; i < 64; ++i) {
  if ((b >> i) & 1) {
   lo ^= a << i;
   if (i) hi ^= a >> (64 - i);
  }
 }
 return {lo, hi};
}
static inline u64 reduce_naive(u64 lo, u64 hi) {
 for (int i = 63; i >= 0; --i) {
  if ((hi >> i) & 1) {
   hi ^= u64(1) << i;
   lo ^= IRRED_LOW << i;
   if (i > 0) hi ^= IRRED_LOW >> (64 - i);
  }
 }
 return lo;
}
static inline u64 mul(u64 a, u64 b) { auto [l, h] = clmul_loop(a, b); return reduce_naive(l, h); }
static inline u64 pow_(u64 a, u64 e) {
 u64 r = 1;
 while (e) { if (e & 1) r = mul(r, a); a = mul(a, a); e >>= 1; }
 return r;
}

int main(int argc, char* argv[]) {
 registerTestlibCmd(argc, argv);
 int T = inf.readInt();
 std::vector<u64> xs(T);
 for (int i = 0; i < T; ++i) xs[i] = (u64) inf.readUnsignedLong();

 for (int i = 0; i < T; ++i) {
  u64 k_user = (u64) ouf.readUnsignedLong();
  u64 v = pow_(G, k_user);
  if (v != xs[i]) {
   quitf(_wa, "case %d: x = %llu, user k = %llu, but g^k = %llu (expected %llu)",
         i, (unsigned long long) xs[i], (unsigned long long) k_user,
         (unsigned long long) v, (unsigned long long) xs[i]);
  }
 }
 // 期待出力の読み飛ばし (内容は確認しないが、Format Error 検知のため)
 for (int i = 0; i < T; ++i) (void) ans.readUnsignedLong();
 quitf(_ok, "OK %d cases", T);
}
