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

// yosupo discrete_logarithm_fixed_mod:
//   p = 999999937 (固定の素数), g (generator), n (queries)
//   各 query で 0 ≤ a < p が与えられ、log_g(a) mod (p-1) を返す
//   (a == 0 → -1、a == 1 → 0)
//
// Library Checker URL:
//   https://judge.yosupo.jp/problem/discrete_logarithm_fixed_mod
