// harness: 各 algos/*.hpp が定義する struct CharPoly::run(N, M) を計測する。
// yosupo "Characteristic Polynomial" 形式の I/O。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {
    cin.tie(0);
    ios::sync_with_stdio(false);
    int N;
    cin >> N;
    vector<vector<u32>> M(N, vector<u32>(N));
    for (auto& row : M) for (auto& x : row) cin >> x;

    constexpr int REPEAT = 1;
    uint64_t best_ns = ~uint64_t(0);
    vector<u32> result;

    for (int rep = 0; rep < REPEAT; ++rep) {
        auto t0 = chrono::steady_clock::now();
        auto r = CharPoly::run(N, M);
        auto t1 = chrono::steady_clock::now();
        result = std::move(r);
        auto ns = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        if (ns < best_ns) best_ns = ns;
    }

    fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
    for (int i = 0; i <= N; ++i) {
        cout << result[i];
        cout << (i == N ? '\n' : ' ');
    }
    return 0;
}
