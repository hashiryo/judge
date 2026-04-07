// judge: https://judge.yosupo.jp/problem/unionfind
#include <bits/stdc++.h>
using namespace std;

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

int main() {
    int n, q;
    scanf("%d %d", &n, &q);
    UnionFind uf(n);
    while (q--) {
        int t, u, v;
        scanf("%d %d %d", &t, &u, &v);
        if (t == 0) uf.unite(u, v);
        else printf("%d\n", uf.same(u, v) ? 1 : 0);
    }
}
