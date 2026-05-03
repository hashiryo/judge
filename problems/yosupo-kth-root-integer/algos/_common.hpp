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

// yosupo Kth Root Integer (https://judge.yosupo.jp/problem/kth_root_integer)
//
// 入力フォーマット:
//   T
//   A_1 k_1
//   A_2 k_2
//   ...
//   A_T k_T
// 出力: floor(A^{1/k}) を各クエリで出力。
//
// 制約: 0 ≤ A ≤ 2^64 - 1, 1 ≤ k ≤ 64
