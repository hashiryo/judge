#pragma once
// algos 共通: typedef とよく使うヘッダ。
// ファイル名が _ 始まりなので detect-changed.py の提出対象列挙からは除外される。
// USE_SIMDE 環境では simde を先に include (float16_t 衝突回避)。
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#endif
#include <bits/stdc++.h>
using namespace std;
