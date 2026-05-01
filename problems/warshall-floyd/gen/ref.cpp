#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>

constexpr std::int32_t WF_INF= 1'000'000'000;

int main() {
 int V;
 std::uint64_t W, seed;
 if(!(std::cin >> V >> W >> seed)) {
  return 0;
 }

 std::vector<std::int32_t> d((std::size_t)V * V, WF_INF);
 constexpr std::uint64_t LCG_MA= 6364136223846793005ULL;
 constexpr std::uint64_t LCG_MB= 1442695040888963407ULL;
 std::uint64_t s= seed;
 for(int i= 0; i < V; ++i) {
  d[(std::size_t)i * V + i]= 0;
  for(int j= 0; j < V; ++j) {
   if(i == j) continue;
   s= s * LCG_MA + LCG_MB;
   d[(std::size_t)i * V + j]= (std::int32_t)((s >> 16) % W) + 1;
  }
 }

 for(int k= 0; k < V; ++k)
  for(int i= 0; i < V; ++i)
   for(int j= 0; j < V; ++j) {
    std::int32_t nd= d[(std::size_t)i * V + k] + d[(std::size_t)k * V + j];
    if(d[(std::size_t)i * V + j] > nd) d[(std::size_t)i * V + j]= nd;
   }

 std::uint64_t acc= 0;
 for(int i= 0; i < V; ++i)
  for(int j= 0; j < V; ++j) acc^= (std::uint64_t)(std::uint32_t)d[(std::size_t)i * V + j];

 std::cout << acc << '\n';
 return 0;
}
