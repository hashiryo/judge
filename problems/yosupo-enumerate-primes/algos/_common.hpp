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

// yosupo Enumerate Primes (https://judge.yosupo.jp/problem/enumerate_primes)
//
// 入力: N A B (N ≤ 5·10^8, 1 ≤ A ≤ π(N), 0 ≤ B < A)
// 出力: 1 行目 "π(N) X"  (X = N 以下の素数のうち index B, B+A, B+2A, ... の個数)
//       2 行目: それら X 個の素数を空白区切りで
