#include <cstdint>
#include <iostream>
int main() {
 std::uint64_t n;
 std::uint64_t mod;
 std::uint64_t state;
 std::uint64_t a;
 std::uint64_t b;
 if(!(std::cin >> n >> mod >> state >> a >> b)) {
  return 0;
 }

 for(std::uint64_t iteration= 0; iteration < n; ++iteration) {
  state= (static_cast<unsigned __int128>(state) * a + b) % mod;
 }

 std::cout << state << '\n';
 return 0;
}
