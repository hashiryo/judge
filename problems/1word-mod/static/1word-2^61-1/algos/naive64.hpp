#pragma once
#include "_common.hpp"
struct DIV {
 constexpr DIV() {}
 constexpr u64 mod(u64 a) const { return a % B; }
};
