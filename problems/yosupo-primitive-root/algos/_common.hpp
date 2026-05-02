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

// Q ≤ 100 個の素数 p (2 ≤ p ≤ 1e18 < 2^62) に対し、原始根 (の 1 つ) を返す。
// 戻り値: 各クエリ p について 1 ≤ g < p を満たす原始根 g。
