#include <cstdint>
#include <iostream>
#include <numeric>
int main() {
 std::uint64_t n;
 std::uint64_t abits;
 std::uint64_t sa;
 std::uint64_t sb;
 if(!(std::cin >> n >> abits >> sa >> sb)) {
  return 0;
 }

 constexpr std::uint64_t LCG_MA= 6364136223846793005ULL;
 constexpr std::uint64_t LCG_MB= 1442695040888963407ULL;
 const std::uint64_t mask= abits >= 64 ? ~std::uint64_t(0) : ((std::uint64_t(1) << abits) - 1);

 std::uint64_t acc= 0;
 for(std::uint64_t i= 0; i < n; ++i) {
  sa= sa * LCG_MA + LCG_MB;
  sb= sb * LCG_MA + LCG_MB;
  std::uint64_t a= sa & mask;
  std::uint64_t b= sb & mask;
  acc^= std::gcd(a, b);
 }

 std::cout << acc << '\n';
 return 0;
}
