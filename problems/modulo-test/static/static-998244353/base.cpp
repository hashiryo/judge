// harness: 各 algos/*.hpp が定義する struct MP を使って計測する
#include <bits/stdc++.h>
using namespace std;
using u32 = unsigned;
using u64 = unsigned long long;
using u128 = __uint128_t;

#include ALGO_HPP   // -DALGO_HPP="\"algos/xxx.hpp\""

signed main() {
    cin.tie(0);
    ios::sync_with_stdio(false);
    constexpr u32 MOD = 998244353;
    u32 n_, state_, a_, b_, c_, d_;
    cin >> n_ >> state_ >> a_ >> b_ >> c_ >> d_;

    // REPEAT は wall time を REPEAT 倍するので TLE と相談して決める。
    // ここは重ループで 1 回でも十分長く動くので REPEAT=1 で十分な精度。
    constexpr int REPEAT = 1;
    uint64_t best_ns = ~uint64_t(0);
    u64 result_out = 0;

    for (int rep = 0; rep < REPEAT; ++rep) {
        constexpr MP mp(MOD);
        auto state = mp.set(state_);
        auto a = mp.set(a_);
        auto b = mp.set(b_);
        auto c = mp.set(c_);
        auto d = mp.set(d_);
        auto t0 = chrono::steady_clock::now();
        for (u32 i = 0; i < n_; ++i) {
            state = mp.mul(mp.plus(mp.mul(state, a), b), mp.plus(mp.mul(state, c), d));
        }
        auto t1 = chrono::steady_clock::now();
        result_out = mp.get(state);
        auto ns = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        if (ns < best_ns) best_ns = ns;
    }

    fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
    cout << result_out << '\n';
    return 0;
}
