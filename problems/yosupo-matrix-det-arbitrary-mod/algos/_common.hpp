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
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// 任意 mod (1 ≤ M ≤ 10^9) 上の N×N 整数行列の行列式 (mod M)。1 ≤ N ≤ 500。
// 入力値 0 ≤ a_{i,j} < M。出力は det(A) mod M (= 0..M-1)。
