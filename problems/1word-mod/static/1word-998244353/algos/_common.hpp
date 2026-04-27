#pragma once
// algos 共通: typedef とよく使うヘッダ + コンパイル時 modulus B。
// ファイル名が _ 始まりなので detect-changed.py の提出対象列挙からは除外される。
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;
constexpr u32 B = 998244353;  // NTT prime, 30 bit
