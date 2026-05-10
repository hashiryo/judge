#pragma once
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#include <simde/x86/clmul.h>
#include <simde/x86/bmi.h>
#else
#include <immintrin.h>
#endif
#include <bits/stdc++.h>
#ifdef __x86_64__
#define PCLMUL [[gnu::target("pclmul")]]
#define VPCLMUL [[gnu::target("vpclmulqdq")]]
#else
#define PCLMUL
#define VPCLMUL
#endif

using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;

// yosupo Convolution over F_{2^64}
//   入力: n m
//         a_0 ... a_{n-1}
//         b_0 ... b_{m-1}
//   出力: c_0 ... c_{n+m-2}
//   ここで c_k = XOR_{i+j=k} a_i · b_j
//   · は F_{2^64} の積 (X^64 + X^4 + X^3 + X + 1 mod 削減多項式)
//
// アルゴリズム:
//   - 素朴: O(nm) で nested loop
//   - additive FFT (a.k.a. nim FFT, characteristic-2 FFT): O((n+m) log(n+m))
//     b[i+1] = b[i]^2 + b[i] で生成する基底列で Frobenius 構造を活用
