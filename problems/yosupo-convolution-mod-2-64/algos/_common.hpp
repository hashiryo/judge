#pragma once
// algos 共通: typedef とよく使うヘッダ。
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// この問題は mod 2^64 (u64 の自然なラップアラウンド)。
// NTT 系は使えず、複素 FFT で u64 を 4 × int16 に分割するのが定番。
