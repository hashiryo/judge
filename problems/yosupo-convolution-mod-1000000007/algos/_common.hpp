#pragma once
// algos 共通: typedef とよく使うヘッダ。
//
// USE_SIMDE 環境 (ARM CI) では simde を <bits/stdc++.h> より先に include する。
// これにより arm_neon.h が <stdfloat> より先にパースされ、`float16_t` が
// arm_neon.h の ::float16_t としてのみ可視になる。後で `using namespace std;` で
// std::float16_t がグローバルに来ても、parse 済みの arm_neon.h 内の参照は影響を
// 受けない (シンボル解決は参照時に行われる)。
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#include <simde/x86/bmi.h>
#endif
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// この問題は mod 10^9+7。NTT-friendly でないので NTT を直接使えない:
// - 複素 FFT (double) + Karatsuba 3-split
// - 多素数 NTT + Garner CRT
// などで対応する。
constexpr u32 MOD = 1'000'000'007;
