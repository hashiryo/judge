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

// mod 998244353 上の N×N 行列の余因子行列 adj(A)。1 ≤ N ≤ 500。
constexpr u32 MOD = 998244353;
