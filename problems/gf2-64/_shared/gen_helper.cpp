// gen_helper.cpp — gf2-64/*/testcases/gen.py の C++ サポーター。
// Python の bit-by-bit GF(2)[x] 演算が遅いので、scalar C++ で書き換える。
// PCLMUL は依存しない (gen はホスト機の機能を仮定したくない、CI でも安定動作したい)。
//
// I/O フォーマット:
//   stdin 1 行目: <op> <T>
//      op = mul | div | pow | sqrt | inv | log_inv (= "g^k = ?" に対する逆 k → g^k)
//   T 行 (op によって列数が変わる):
//      mul:     a b   →  a*b
//      div:     a b   →  a/b   (b!=0 前提)
//      pow:     a e   →  a^e
//      sqrt:    a     →  sqrt(a) = a^(2^63)
//      inv:     a     →  a^(-1) (a!=0 前提)
//      log_inv: k     →  g^k   (g = 2 固定)
//
// stdout に T 行で結果。
//
// Build: c++ -std=c++17 -O3 -o /tmp/gf2_64_gen_helper gen_helper.cpp
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

using u64 = unsigned long long;
constexpr u64 IRRED_LOW = 0x1Bull;
constexpr u64 G = 2;

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

static inline u64 mul(u64 a, u64 b) {
 auto [lo, hi] = clmul_loop(a, b);
 return reduce_naive(lo, hi);
}
static inline u64 pow_(u64 a, u64 e) {
 u64 r = 1;
 while (e) { if (e & 1) r = mul(r, a); a = mul(a, a); e >>= 1; }
 return r;
}
static inline u64 inv(u64 a) { return pow_(a, ~u64(1)); } // a^(2^64 - 2)
static inline u64 sqrt_(u64 a) { return pow_(a, u64(1) << 63); }

int main(int argc, char* argv[]) {
 (void) argc; (void) argv;
 // 高速 I/O 用バッファ (T が大きい場合に効く)。
 // setvbuf は最初の I/O 前に呼ぶのが必須なので main 冒頭で実行。
 static char in_buf[1 << 20], out_buf[1 << 20];
 setvbuf(stdin, in_buf, _IOFBF, sizeof in_buf);
 setvbuf(stdout, out_buf, _IOFBF, sizeof out_buf);

 char op[32];
 int T = 0;
 if (scanf("%31s %d", op, &T) != 2) return 1;
 std::string opstr(op);

 if (opstr == "mul" || opstr == "div" || opstr == "pow") {
  for (int i = 0; i < T; ++i) {
   u64 a, b;
   if (scanf("%llu %llu", &a, &b) != 2) return 1;
   u64 r;
   if (opstr == "mul") r = mul(a, b);
   else if (opstr == "div") r = mul(a, inv(b));
   else /* pow */         r = pow_(a, b);
   printf("%llu\n", r);
  }
 } else if (opstr == "sqrt" || opstr == "inv" || opstr == "log_inv") {
  for (int i = 0; i < T; ++i) {
   u64 a;
   if (scanf("%llu", &a) != 1) return 1;
   u64 r;
   if (opstr == "sqrt")         r = sqrt_(a);
   else if (opstr == "inv")     r = inv(a);
   else /* log_inv: k → g^k */  r = pow_(G, a);
   printf("%llu\n", r);
  }
 } else {
  fprintf(stderr, "unknown op: %s\n", op);
  return 1;
 }
 return 0;
}
