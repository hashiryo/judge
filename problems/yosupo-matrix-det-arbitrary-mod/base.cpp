// harness: 各 algos/*.hpp が定義する struct Det::run(N, mod, A) を計測する。
// yosupo の "Matrix Determinant (arbitrary mod)" 形式の I/O。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {
    cin.tie(0);
    ios::sync_with_stdio(false);
    int N;
    u32 MOD;
    cin >> N >> MOD;
    vector<vector<u32>> a(N, vector<u32>(N));
    for (auto& row : a)
        for (auto& x : row) cin >> x;

    constexpr int REPEAT = 1;
    uint64_t best_ns = ~uint64_t(0);
    u32 result = 0;

    for (int rep = 0; rep < REPEAT; ++rep) {
        auto t0 = chrono::steady_clock::now();
        u32 r = Det::run(N, MOD, a);
        auto t1 = chrono::steady_clock::now();
        result = r;
        auto ns = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        if (ns < best_ns) best_ns = ns;
    }

    fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
    cout << result << '\n';
    return 0;
}
