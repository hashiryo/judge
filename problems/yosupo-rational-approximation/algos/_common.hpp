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

// yosupo Rational Approximation (https://judge.yosupo.jp/problem/rational_approximation)
//
// 各クエリで (n, x, y) が与えられ、a/b ≤ x/y ≤ c/d で 1 ≤ a,b,c,d ≤ n、
// かつ |a/b - x/y|, |c/d - x/y| が最小となる (a, b, c, d) を出力。
// (Stern-Brocot 木上の二分探索で求まる)
