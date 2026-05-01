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

// mod 998244353 上の N×M 行列 A, 長さ N の b に対し Ax=b の解空間を返す。
// 1 ≤ N, M ≤ 500。
constexpr u32 MOD = 998244353;

// 解の表現:
//   ok=true なら (sol: 長さ M の特解, basis: R 行 M 列の解空間基底)
//   ok=false なら sol/basis は空 (解なし)
struct SolveResult {
 bool ok;
 vector<u32> sol;
 vector<vector<u32>> basis;
};
