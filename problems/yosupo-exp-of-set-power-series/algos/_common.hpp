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

// mod 998244353 上の set power series exp。0 ≤ N ≤ 20、長さ 2^N、b_0 = 0。
//   exp(s) を mod (x_0^2, ..., x_{N-1}^2) で求める。
constexpr u32 MOD = 998244353;
