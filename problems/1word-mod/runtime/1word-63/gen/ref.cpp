#include <cstdint>
#include <iostream>
int main() {
 std::uint64_t n, b_in, sa, ma, mb;
 if(!(std::cin >> n >> b_in >> sa >> ma >> mb)) {
  return 0;
 }
 constexpr std::uint64_t MASK= (std::uint64_t(1) << 63) - 1;
 const std::uint64_t b= b_in;
 std::uint64_t acc= 0;
 for(std::uint64_t i= 0; i < n; ++i) {
  sa= sa * ma + mb;
  std::uint64_t a= sa & MASK;
  acc^= a % b;
 }
 std::cout << acc << '\n';
 return 0;
}
