#pragma once
// algos 共通: typedef とよく使うヘッダ。
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// yosupo convolution_mod の固定 mod。NTT 素数 (998244353 = 119 * 2^23 + 1, ω=3)。
constexpr u32 MOD = 998244353;
