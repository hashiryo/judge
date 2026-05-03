#include <cstdint>
#include <iostream>
int main() {
 std::uint64_t n;
 std::uint64_t mod;
 std::uint64_t bbits;
 std::uint64_t sa;
 std::uint64_t sb;
 if(!(std::cin >> n >> mod >> bbits >> sa >> sb)) {
  return 0;
 }

 constexpr std::uint64_t LCG_MA = 6364136223846793005ULL;
 constexpr std::uint64_t LCG_MB = 1442695040888963407ULL;
 const std::uint32_t b_mask = bbits >= 32 ? ~std::uint32_t(0) : ((std::uint32_t(1) << bbits) - 1);

 std::uint32_t acc = 0;
 for(std::uint64_t i = 0; i < n; ++i) {
  sa = sa * LCG_MA + LCG_MB;
  sb = sb * LCG_MA + LCG_MB;
  std::uint32_t a = (std::uint32_t)(sa % mod);
  std::uint32_t e = (std::uint32_t)sb & b_mask;
  std::uint64_t r = 1 % mod;
  std::uint64_t base = a;
  while(e) {
   if(e & 1) r = r * base % mod;
   base = base * base % mod;
   e >>= 1;
  }
  acc ^= (std::uint32_t)r;
 }

 std::cout << acc << '\n';
 return 0;
}
