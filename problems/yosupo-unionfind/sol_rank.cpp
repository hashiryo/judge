#include <bits/stdc++.h>
using namespace std;

struct UnionFind {
    vector<int> par, rank_;
    UnionFind(int n) : par(n), rank_(n, 0) { iota(par.begin(), par.end(), 0); }
    int find(int x) { return par[x] == x ? x : par[x] = find(par[x]); }
    bool unite(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) return false;
        if (rank_[x] < rank_[y]) swap(x, y);
        par[y] = x;
        if (rank_[x] == rank_[y]) rank_[x]++;
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
