#pragma once
// USE_SIMDE 環境では simde を先に include (float16_t 衝突回避)。
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

// この問題は GF(2) (mod 2) 上の N×M {0,1} 行列のランク。1 ≤ N, M ≤ 4096。
// 結果は 0 以上 min(N, M) 以下の整数。
