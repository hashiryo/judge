// harness: 各 algos/*.hpp が定義する struct Solver::run(N, M, A, b) を計測する。
// yosupo "System of Linear Equations" 形式の I/O。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive.hpp"
#endif
#include ALGO_HPP

signed main() {
    cin.tie(0);
    ios::sync_with_stdio(false);
    int N, M;
    cin >> N >> M;
    vector<vector<u32>> A(N, vector<u32>(M));
    for (auto& row : A) for (auto& x : row) cin >> x;
    vector<u32> b(N);
    for (auto& x : b) cin >> x;

    constexpr int REPEAT = 1;
    uint64_t best_ns = ~uint64_t(0);
    SolveResult result;

    for (int rep = 0; rep < REPEAT; ++rep) {
        auto t0 = chrono::steady_clock::now();
        auto r = Solver::run(N, M, A, b);
        auto t1 = chrono::steady_clock::now();
        result = std::move(r);
        auto ns = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        if (ns < best_ns) best_ns = ns;
    }

    fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
    if (!result.ok) {
        cout << -1 << '\n';
    } else {
        cout << result.basis.size() << '\n';
        for (size_t j = 0; j < result.sol.size(); ++j) {
            cout << result.sol[j];
            cout << (j + 1 == result.sol.size() ? '\n' : ' ');
        }
        for (auto const& row : result.basis) {
            for (size_t j = 0; j < row.size(); ++j) {
                cout << row[j];
                cout << (j + 1 == row.size() ? '\n' : ' ');
            }
        }
    }
    return 0;
}
