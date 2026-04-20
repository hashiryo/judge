#include <cstdint>
#include <iostream>
namespace {

constexpr std::uint64_t MOD= 1000000007;

}  // namespace
int main() {
 std::uint64_t n;
 std::uint64_t state;
 std::uint64_t a;
 std::uint64_t b;
 std::uint64_t c;
 std::uint64_t d;
 if(!(std::cin >> n >> state >> a >> b >> c >> d)) {
  return 0;
 }

 for(std::uint64_t iteration= 0; iteration < n; ++iteration) {
  const std::uint64_t t1= (static_cast<unsigned __int128>(state) * a + b) % MOD;
  const std::uint64_t t2= (static_cast<unsigned __int128>(state) * c + d) % MOD;
  state= static_cast<unsigned __int128>(t1) * t2 % MOD;
 }

 std::cout << state << '\n';
 return 0;
}
