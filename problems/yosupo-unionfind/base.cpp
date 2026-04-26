// harness: 各 algos/*.hpp が定義する struct UnionFind を使って計測する
#include "algos/_common.hpp"

// CI では -DALGO_HPP="\"algos/xxx.hpp\"" で上書きされる。
// ここでのデフォルトは IDE で base.cpp 単独表示時に補完を効かせるため。
#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {
    int n, q;
    if (scanf("%d %d", &n, &q) != 2) return 1;
    vector<array<int, 3>> qs(q);
    for (auto& e : qs) scanf("%d %d %d", &e[0], &e[1], &e[2]);

    // クエリ列を先読みしてあるので、計測対象は純粋な UnionFind 演算のみ。
    UnionFind uf(n);
    string out;
    out.reserve(size_t(q) * 2);
    auto t0 = chrono::steady_clock::now();
    for (auto& e : qs) {
        if (e[0] == 0) uf.unite(e[1], e[2]);
        else out += uf.same(e[1], e[2]) ? "1\n" : "0\n";
    }
    auto t1 = chrono::steady_clock::now();
    auto ns = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();

    fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)ns);
    fputs(out.c_str(), stdout);
    return 0;
}
