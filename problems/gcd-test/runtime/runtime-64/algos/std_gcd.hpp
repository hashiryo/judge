#pragma once
#include "_common.hpp"
#include <numeric>
// libstdc++ の std::gcd。中身は __detail::__gcd の Euclidean ベース実装。
struct G {
 static inline u64 gcd(u64 a, u64 b) { return std::gcd(a, b); }
};
