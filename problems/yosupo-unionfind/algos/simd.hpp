#pragma once
#include "_common.hpp"
// SIMD の動作確認用。Union-Find 自体は SIMD の恩恵が薄いが、
// ifdef の切り替えとコンパイルが通ることの確認として使う。
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#else
#include <immintrin.h>
#endif

struct UnionFind {
    vector<int> par;
    UnionFind(int n) : par(n, -1) {
        // SIMD コンパイル確認: AVX2 の命令が使えるか
        __m256i a = _mm256_set1_epi32(42);
        __m256i b = _mm256_set1_epi32(1);
        __m256i c = _mm256_add_epi32(a, b);
        (void)c;
    }
    int find(int x) { return par[x] < 0 ? x : par[x] = find(par[x]); }
    bool unite(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) return false;
        if (par[x] > par[y]) swap(x, y);
        par[x] += par[y];
        par[y] = x;
        return true;
    }
    bool same(int x, int y) { return find(x) == find(y); }
};
