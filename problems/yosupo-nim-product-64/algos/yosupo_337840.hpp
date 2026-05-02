#pragma once
// =============================================================================
// Source: yosupo "Nim Product (F_{2^64})" 提出 337840.
//   https://judge.yosupo.jp/submission/337840
//   方式:
//     - 線形同型 F_{2^64} ≅ GF(2)[x]/(x^64 + x^4 + x^3 + x + 1) を介する。
//       nim 表現 → poly 表現は INV_COL[64] (基底変換行列) のバイト分割 lookup。
//       poly 表現 → nim 表現は BASIS_COL[64] の同様の lookup。
//     - poly × poly は PCLMUL (`_mm_clmulepi64_si128`) で 1 命令。
//     - 128-bit 結果を上記原始多項式で reduce (上位 4 bit 用テーブルを併用)。
//   結果: 1 nim 積 = 8 byte lookup × 16 + 1 PCLMUL + 数 cycle reduce。
// 抽出方針:
//   - 元の I/O (FastIO の mmap stdin / writebuf stdout) と big_alloc は省略
//     (base.cpp 側の cin/cout で揃える)
//   - BASIS_COL / INV_COL / byte tables / clmul / reduce_mod / nim 全変換 +
//     nim_mul を namespace ラップでそのまま流用
//   - run() で 4-pass (nim→poly, clmul+reduce, poly→nim, write) に分けて
//     OoO に乗りやすくする (cp-algo の checkpoint 構造と同じ)
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt,pclmul")
#endif
#include "_common.hpp"
#ifdef USE_SIMDE
// arm でも SIMDe が x86 intrinsic を提供してくれる (内部で NEON / scalar fallback)。
#include <simde/x86/sse2.h>
#include <simde/x86/clmul.h>
#else
#include <immintrin.h>
#endif

namespace yosupo_337840 {
using std::array;
using std::uint64_t;
using std::uint8_t;

constexpr auto BASIS_COL = array{
 0x0000000000000001ull, 0x5211145c804b6109ull, 0x7c8bc2cad259879full, 0x565854b4c60c1e0bull,
 0x4068acf7104c20c3ull, 0x662d2bd0f2739155ull, 0x7a90c83701fa8323ull, 0x21cfa750247e8755ull,
 0x67d1044e545abf47ull, 0x4d9d3b5a8568f839ull, 0x567a9d7331b6b3c6ull, 0x1ca54bfdd6d1ae59ull,
 0x454fa483275db25cull, 0x6766df6fec4e9d44ull, 0x35cb621cec1fe7f9ull, 0x4c606d3e52faf263ull,
 0x57640dc825a57954ull, 0x7aca87838b7f6315ull, 0x6d53c884ebf2b0edull, 0x3721d998bb50164bull,
 0x7aa7c62fd6cd53abull, 0x47cbb2c51f7c040full, 0x132063b7f5e42489ull, 0x0c1b36c8b2993f8aull,
 0x60119ecff680497aull, 0x5175da444cc11791ull, 0x5792ff4554765b09ull, 0x0c9fdb8a01334e82ull,
 0x2be0a763a68a4725ull, 0x3c2dc8260ad051f6ull, 0x6c4c9fed8816bb9cull, 0x630062753ffaf766ull,
 0x7b37d31b5d519225ull, 0x2364f7f79705691cull, 0x453eb8a83e2fec71ull, 0x7c0121b37e828666ull,
 0x59190d3250e66011ull, 0x103207f9dda18caeull, 0x28233dce01c69b76ull, 0x4fa519899227a5e7ull,
 0x4567ba46ee7bc6cdull, 0x0a284773d021afd5ull, 0x63894079bbe3a824ull, 0x11013c7fdfaaa5c2ull,
 0x1aa984f18574f3b0ull, 0x0cbaba126fd0c4dbull, 0x0b8797719e6dc725ull, 0x4a2845680aefaa72ull,
 0x536d2535f6934e15ull, 0x01db7a57effcd689ull, 0x7e1ed0ad01e2a5adull, 0x0aedc9b3cee826f6ull,
 0x7ba716eccf9f68e1ull, 0x5d5e23bc0f3dc38full, 0x0b5f2a3b88674d83ull, 0x2de9bafc2f00f8d4ull,
 0x3b56712ad419c7e0ull, 0x3ab4be8c30c19253ull, 0x2708522ffaa654b0ull, 0x2b8bca57bf643598ull,
 0x588825d1a5fa8e1cull, 0x86adf8bf4d45962full, 0x51b4c15d8719dd73ull, 0xe4a2b3b59783d0aaull
};

constexpr auto INV_COL = array{
 0x0000000000000001ull, 0x19c9369f278adc02ull, 0xa181e7d66f5ff795ull, 0x5db84357ce785d09ull,
 0xa0bae2f9d2430cc8ull, 0xb7ea5a9705b771c0ull, 0xba4f3cd82801769dull, 0x4886cde01b8241d0ull,
 0x0a6f43f2aaf612edull, 0xebd0142f98030a32ull, 0xa81f89cda43f3792ull, 0xe99aec6b66ccb814ull,
 0xa69d1ff025fc2f82ull, 0x48a81132d25db068ull, 0x4a900f9dcaa9644full, 0xe5ce4ea88259972aull,
 0xf7094c336029f04cull, 0xe191dde287bc9c6bull, 0xaacaff12bff239b8ull, 0x49bc5212be1bc1caull,
 0xfe57defb454446cfull, 0xa1dffcf944bdf6a7ull, 0xb9f1bdb5cee941eeull, 0x12e5e889275c22deull,
 0x5bcb6b117b77eeedull, 0x03eb1ab59d05ae4bull, 0x02a25d7076ddd386ull, 0x53164a606c612245ull,
 0xebb33f5822f66059ull, 0xe9be765f5747b93eull, 0x552a78df373a354full, 0xbcf5ac65f31fb8bfull,
 0xe411e728becdc77bull, 0xf35c26d7b57cdca6ull, 0x4499da83de4ca5f7ull, 0x40ab25bdca4ae226ull,
 0xee004b6f1dff7218ull, 0x0d122da9821c5b41ull, 0x51fbfcb058120efeull, 0xa148b1fa84905b22ull,
 0xbb8ed3e647604d8dull, 0xe2d93fef2472776full, 0x4c17a2541a10e6b5ull, 0x1d879e08903708e7ull,
 0x0fbe7d0d1934da90ull, 0x5bf977d9c6f61d30ull, 0x06832fc918260412ull, 0x0fe22e843ebf73e3ull,
 0x4d7ef4e4fa28d60dull, 0x402250d979afbed5ull, 0x067902b8c8ca2d4full, 0xf38d113fe1d6bb16ull,
 0x414f0248b02b5b7dull, 0xf041922915824ce9ull, 0x11a72fb5e30c93d9ull, 0x12e54f4d63102aeeull,
 0xbc46ac14b3141c6cull, 0x1f172b3c16c645bbull, 0x584b492ed4e8fa6cull, 0x00a852e9a32cc133ull,
 0xa180861bce00a45eull, 0xa194b6bcb4645fb9ull, 0x4509002ad808a4fbull, 0xc5172a0055602f69ull
};

template <const auto& COLS>
constexpr auto make_byte_tables() {
 array<array<uint64_t, 1 << 8>, 8> T{};
 for (int pos = 0; pos < 8; ++pos) {
  for (int col = 0; col < 8; ++col) {
   for (int mask = 0; mask < (1 << col); ++mask) {
    T[pos][mask | (1 << col)] = T[pos][mask] ^ COLS[pos * 8 + col];
   }
  }
 }
 return T;
}

static constexpr auto INV_BYTE = make_byte_tables<INV_COL>();
static constexpr auto BASIS_BYTE = make_byte_tables<BASIS_COL>();

[[gnu::always_inline]]
inline uint64_t nim_to_poly(uint64_t x) {
 auto xb = std::bit_cast<array<uint8_t, 8>>(x);
 return INV_BYTE[0][xb[0]] ^ INV_BYTE[1][xb[1]]
      ^ INV_BYTE[2][xb[2]] ^ INV_BYTE[3][xb[3]]
      ^ INV_BYTE[4][xb[4]] ^ INV_BYTE[5][xb[5]]
      ^ INV_BYTE[6][xb[6]] ^ INV_BYTE[7][xb[7]];
}

[[gnu::always_inline]]
inline uint64_t poly_to_nim(uint64_t c) {
 auto cb = std::bit_cast<array<uint8_t, 8>>(c);
 return BASIS_BYTE[0][cb[0]] ^ BASIS_BYTE[1][cb[1]]
      ^ BASIS_BYTE[2][cb[2]] ^ BASIS_BYTE[3][cb[3]]
      ^ BASIS_BYTE[4][cb[4]] ^ BASIS_BYTE[5][cb[5]]
      ^ BASIS_BYTE[6][cb[6]] ^ BASIS_BYTE[7][cb[7]];
}

// Carryless multiply over GF(2)
[[gnu::always_inline]]
inline __m128i clmul(uint64_t a, uint64_t b) {
 // SIMDe は _mm_clmulepi64_si128 を関数マクロにしているので、引数の braced-init を
 // そのまま渡せない (マクロ展開でカンマが余分な引数として解釈される)。一時変数経由。
 __m128i av{(long long)a, 0};
 __m128i bv{(long long)b, 0};
 return _mm_clmulepi64_si128(av, bv, 0);
}

// Reduce modulo x^64 + x^4 + x^3 + x + 1
[[gnu::always_inline]]
inline uint64_t reduce_mod(__m128i v) {
 v[0] ^= v[1] ^ (v[1] << 1) ^ (v[1] << 3) ^ (v[1] << 4);
 static constexpr auto RED_OVER = [] {
  array<uint64_t, 16> red{};
  for (int q = 0; q < 16; ++q) {
   uint64_t o = q ^ (q >> 1) ^ (q >> 3);
   red[q] = o ^ (o << 1) ^ (o << 3) ^ (o << 4);
  }
  return red;
 }();
 return v[0] ^ RED_OVER[v[1] >> 60];
}

[[gnu::always_inline]]
inline uint64_t nim_mul(uint64_t a, uint64_t b) {
 uint64_t pa = nim_to_poly(a);
 uint64_t pb = nim_to_poly(b);
 auto prod = clmul(pa, pb);
 uint64_t red = reduce_mod(prod);
 return poly_to_nim(red);
}

} // namespace yosupo_337840

struct NimProduct {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using namespace yosupo_337840;
  const size_t T = as.size();
  vector<u64> ans(T);
  // 元提出は 4-pass (nim→poly, clmul+reduce, poly→nim, write) で
  // OoO 並列化を狙う構造。同じ形に保つ。
  vector<u64> A(T), B(T);
  for (size_t i = 0; i < T; ++i) A[i] = nim_to_poly(as[i]);
  for (size_t i = 0; i < T; ++i) B[i] = nim_to_poly(bs[i]);
  for (size_t i = 0; i < T; ++i) A[i] = reduce_mod(clmul(A[i], B[i]));
  for (size_t i = 0; i < T; ++i) ans[i] = poly_to_nim(A[i]);
  return ans;
 }
};
