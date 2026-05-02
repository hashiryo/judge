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

// T 個の (a, b) ∈ F_{2^64}^2 に対し nim 積 a ⊗ b を計算する。
// T ≤ 1e6。出力一意 (= checker 不要、diff で判定可能)。
