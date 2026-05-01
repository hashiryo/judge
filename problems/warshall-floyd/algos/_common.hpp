#pragma once
// algos 共通: typedef とよく使うヘッダ。
// USE_SIMDE 環境では simde を先に include (float16_t 衝突回避)。
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#endif
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i32 = int;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// distance 行列の表現: V*Vp の i32 (Vp = V を 8 の倍数に切り上げた値)。
// 未到達は INF。INF + 任意の有効距離 が i32 に収まるよう INF=10^9 を採用 (2*INF < 2^31)。
constexpr i32 WF_INF = 1'000'000'000;
