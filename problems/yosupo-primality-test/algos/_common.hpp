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

// 64-bit 整数 N (1 ≤ N ≤ 1e18) に対する素数判定。Q ≤ 1e5。
// 戻り値: true/false の vector を返す (base.cpp が "Yes"/"No" を出力)。
