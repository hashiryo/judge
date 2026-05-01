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

// mod 998244353 上の polynomial composite set power series。
//   入力: f (長さ M, 一変数多項式), b (長さ 2^N, set power series)
//   出力: c = f(b) mod (x_0^2, ..., x_{N-1}^2)
// 0 ≤ M ≤ 1e5, 0 ≤ N ≤ 20。
constexpr u32 MOD = 998244353;
