#include <cstdint>
#include <iostream>
int main() {
 std::uint64_t n, sa, ma, mb;
 if(!(std::cin >> n >> sa >> ma >> mb)) {
  return 0;
 }
 constexpr std::uint32_t MASK= (std::uint32_t(1) << 31) - 1;
 constexpr std::uint32_t B= 998244353;
 std::uint32_t acc= 0;
 for(std::uint64_t i= 0; i < n; ++i) {
  sa= sa * ma + mb;
  std::uint32_t a= std::uint32_t(sa) & MASK;
  acc^= a % B;
 }
 std::cout << acc << '\n';
 return 0;
}
