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

// mod 998244353 上の N×N 行列の特性多項式 p(x) = det(xI - M) を求める。
// 0 ≤ N ≤ 500。出力は p_0, p_1, ..., p_N の N+1 係数 (p_N = 1, モニック)。
constexpr u32 MOD = 998244353;
