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
 if(!(std::cin >> n >> state >> a >> b)) {
  return 0;
 }

 for(std::uint64_t iteration= 0; iteration < n; ++iteration) {
  state= (static_cast<unsigned __int128>(state) * a + b) % MOD;
 }

 std::cout << state << '\n';
 return 0;
}
