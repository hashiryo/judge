#pragma once
// algos 共通: typedef とよく使うヘッダ。
//
// USE_SIMDE 環境 (ARM CI) では simde を <bits/stdc++.h> より先に include する。
// arm_neon.h の ::float16_t と <stdfloat> の std::float16_t の衝突を防ぐ。
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

// この問題は GF(2) (mod 2) 上の N×N {0,1} 行列の行列式。1 ≤ N ≤ 4096。
// 結果は 0 または 1。
