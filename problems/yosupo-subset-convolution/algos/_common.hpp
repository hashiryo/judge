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

// mod 998244353 上の subset convolution。0 ≤ N ≤ 20、長さ 2^N。
//   c_k = sum_{i AND j = 0, i OR j = k} a_i b_j
constexpr u32 MOD = 998244353;
