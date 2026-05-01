#pragma once
// USE_SIMDE 環境では simde を先に include (float16_t 衝突回避)。
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

// 任意 mod 998244353 上の N×N 整数行列の行列式 (mod 998244353)。N ≤ 500。
constexpr u32 MOD = 998244353;
