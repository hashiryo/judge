#pragma once
#include "_common.hpp"
struct DIV {
 u64 b;
 constexpr DIV(u64 b): b(b) {}
 constexpr u64 mod(u64 a) const { return a % b; }
};
