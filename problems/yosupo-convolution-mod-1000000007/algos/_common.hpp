#pragma once
// algos 共通: typedef とよく使うヘッダ。
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// この問題は mod 10^9+7。NTT-friendly でないので NTT を直接使えない:
// - 複素 FFT (double) + Karatsuba 3-split
// - 多素数 NTT + Garner CRT
// などで対応する。
constexpr u32 MOD = 1'000'000'007;
