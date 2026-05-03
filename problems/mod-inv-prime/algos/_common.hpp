#pragma once
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#endif
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i32 = int32_t;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// 自前実験: p = 998244353 (NTT 素数) 固定で、n 個の a に対する a^{-1} mod p を求める。
//   入力: p n / a_1 / ... / a_n
//   出力: a_1^{-1}, ..., a_n^{-1}
//   a == 0 のとき: -1 を出力 (sentinel)
//
// 関連:
//   - Maspy 記事: https://maspypy.com/o1-mod-inv-mod-pow
//   - モジュラ逆数全般: https://zenn.dev/mizar/articles/ea1aaa320a4b9d
