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

// mod 998244353 上の set power series power projection。
//   入力: a (長さ 2^N), w (長さ 2^N), M
//   出力: ans_m = sum_i w_i * [a^m]_i  (m = 0..M-1)
// 0 ≤ N ≤ 20, 0 ≤ M ≤ 1e5。
constexpr u32 MOD = 998244353;
