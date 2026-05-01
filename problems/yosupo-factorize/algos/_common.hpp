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

// Q ≤ 100 個の数 a (1 ≤ a ≤ 1e18) を素因数分解する。
// 戻り値: 各 i に対し ascending な因数列 (重複あり)。
