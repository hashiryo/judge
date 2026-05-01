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

// mod 998244353 上の N×N 整数行列の逆行列。1 ≤ N ≤ 500。
constexpr u32 MOD = 998244353;

// 戻り値: ok=false なら逆行列なし、ok=true なら mat[N][N] が逆行列 (mod MOD)。
struct InverseResult {
 bool ok;
 vector<vector<u32>> mat;
};
