#pragma once
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#endif
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// 既約多項式 P(x) = x^64 + x^4 + x^3 + x + 1 (= 0x1B + x^64) を fix。
// 上位 4 bit を含まない 64-bit 表現での "x^64 mod P" の下位 64 bit:
constexpr uint64_t IRRED_LOW = 0x1Bu;
