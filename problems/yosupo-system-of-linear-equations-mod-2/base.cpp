// harness: 各 algos/*.hpp が定義する struct Solve::run(N, M, A, b) を計測する。
// yosupo の "System of Linear Equations (mod 2)" 形式の I/O。
// run() が空 vector を返したら解なし → "-1" を出力。
// それ以外は最初に R (= 解空間次元 = 戻り vector のサイズ-1) を出力、続いて R+1 行。
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
    vector<string> a(N);
    for (auto& s : a) cin >> s;
    // b は N 個の '0'/'1' を連結した文字列 (yosupo 形式)。
    string b_str;
    cin >> b_str;
    vector<int> b(N);
    for (int i = 0; i < N; ++i) b[i] = b_str[i] - '0';

    constexpr int REPEAT = 1;
    uint64_t best_ns = ~uint64_t(0);
    vector<string> result;

    for (int rep = 0; rep < REPEAT; ++rep) {
        auto t0 = chrono::steady_clock::now();
        auto r = Solve::run(N, M, a, b);
        auto t1 = chrono::steady_clock::now();
        result = std::move(r);
        auto ns = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        if (ns < best_ns) best_ns = ns;
    }

    fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
    if (result.empty()) cout << "-1\n";
    else {
        cout << (result.size() - 1) << '\n';
        for (const auto& s : result) cout << s << '\n';
    }
    return 0;
}
