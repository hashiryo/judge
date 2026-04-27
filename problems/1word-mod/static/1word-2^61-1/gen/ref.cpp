#include <cstdint>
#include <iostream>
int main() {
 std::uint64_t n, sa, ma, mb;
 if(!(std::cin >> n >> sa >> ma >> mb)) {
  return 0;
 }
 constexpr std::uint64_t MASK= (std::uint64_t(1) << 63) - 1;
 constexpr std::uint64_t B= (std::uint64_t(1) << 61) - 1;
 std::uint64_t acc= 0;
 for(std::uint64_t i= 0; i < n; ++i) {
  sa= sa * ma + mb;
  std::uint64_t a= sa & MASK;
  acc^= a % B;
 }
 std::cout << acc << '\n';
 return 0;
}
