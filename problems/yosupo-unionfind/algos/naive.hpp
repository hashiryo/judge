#pragma once
#include "_common.hpp"
// par[root] = -size, par[child] = parent. 経路圧縮 + サイズ union by size。
struct UnionFind {
    vector<int> par;
    UnionFind(int n) : par(n, -1) {}
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
