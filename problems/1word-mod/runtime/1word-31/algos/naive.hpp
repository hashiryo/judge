#pragma once
#include "_common.hpp"
struct DIV {
 u32 b;
 constexpr DIV(u32 b): b(b) {}
 constexpr u32 mod(u32 a) const { return a % b; }
};
