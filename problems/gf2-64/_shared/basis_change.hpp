#pragma once
// =============================================================================
// poly 基底 (GF(2)[x] / (x^64+x^4+x^3+x+1)) ↔ nim 自然表現 の基底変換。
// yosupo_337840.hpp で使われているテーブル (BASIS_COL / INV_COL) と同等。
// 既約多項式 P(x) = x^64+x^4+x^3+x+1 なので、yosupo の表をそのまま流用できる。
//
// 用途: 塔基底アルゴリズム (tower.hpp) で poly 入力を nim 表現に持ち込み、
//       mylib の Nimber を呼ぶ → 戻すのに使う。
// =============================================================================
#include <array>
#include <bit>
#include <cstdint>

namespace gf2_64_basis {
using u64 = unsigned long long;
using u8 = unsigned char;

// poly → nim の基底変換用 64 個の "列ベクトル" (各 u64 は 64-bit GF(2) ベクトル)
constexpr std::array<u64, 64> BASIS_COL = {
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

// nim → poly の基底変換用 (BASIS_COL の逆写像)
constexpr std::array<u64, 64> INV_COL = {
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

// 8 個の byte → u64 ルックアップテーブルにコンパイル時展開
template <const auto& COLS>
constexpr auto make_byte_tables() {
 std::array<std::array<u64, 256>, 8> T{};
 for (int pos = 0; pos < 8; ++pos)
  for (int col = 0; col < 8; ++col)
   for (int mask = 0; mask < (1 << col); ++mask)
    T[pos][mask | (1 << col)] = T[pos][mask] ^ COLS[pos * 8 + col];
 return T;
}
inline constexpr auto BASIS_BYTE = make_byte_tables<BASIS_COL>();
inline constexpr auto INV_BYTE = make_byte_tables<INV_COL>();

// poly 基底 u64 → nim 自然表現 u64
[[gnu::always_inline]] inline u64 poly_to_nim(u64 x) {
 auto xb = std::bit_cast<std::array<u8, 8>>(x);
 return BASIS_BYTE[0][xb[0]] ^ BASIS_BYTE[1][xb[1]]
      ^ BASIS_BYTE[2][xb[2]] ^ BASIS_BYTE[3][xb[3]]
      ^ BASIS_BYTE[4][xb[4]] ^ BASIS_BYTE[5][xb[5]]
      ^ BASIS_BYTE[6][xb[6]] ^ BASIS_BYTE[7][xb[7]];
}

// nim 自然表現 u64 → poly 基底 u64
[[gnu::always_inline]] inline u64 nim_to_poly(u64 x) {
 auto xb = std::bit_cast<std::array<u8, 8>>(x);
 return INV_BYTE[0][xb[0]] ^ INV_BYTE[1][xb[1]]
      ^ INV_BYTE[2][xb[2]] ^ INV_BYTE[3][xb[3]]
      ^ INV_BYTE[4][xb[4]] ^ INV_BYTE[5][xb[5]]
      ^ INV_BYTE[6][xb[6]] ^ INV_BYTE[7][xb[7]];
}

} // namespace gf2_64_basis
