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

constexpr uint64_t IRRED_LOW = 0x1Bu;
constexpr uint64_t LOG_GENERATOR = 2;  // 原始元 (P(x) = x^64+x^4+x^3+x+1 の下で)
constexpr uint64_t GROUP_ORDER = 0xFFFFFFFFFFFFFFFFull;  // = 2^64 - 1
