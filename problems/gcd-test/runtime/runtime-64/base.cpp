// harness: 各 algos/*.hpp が定義する struct G の static G::gcd(u64, u64) を計測する。
#include "algos/_common.hpp"

// CI では -DALGO_HPP="\"algos/xxx.hpp\"" で上書きされる。
// ここでのデフォルトは IDE で base.cpp 単独表示時に補完を効かせるため。
#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive_mod.hpp"
#endif
#include ALGO_HPP

signed main() {
    cin.tie(0);
    ios::sync_with_stdio(false);
    u64 n_, abits_, sa_, sb_;
    cin >> n_ >> abits_ >> sa_ >> sb_;

    // (a_i, b_i) を計測前に LCG で生成しておく (生成コストを algo time から除外)。
    constexpr u64 LCG_MA = 6364136223846793005ULL;
    constexpr u64 LCG_MB = 1442695040888963407ULL;
    const u64 mask = abits_ >= 64 ? ~u64(0) : ((u64(1) << abits_) - 1);
    vector<u64> as(n_), bs(n_);
    {
        u64 sa = sa_, sb = sb_;
        for (u64 i = 0; i < n_; ++i) {
            sa = sa * LCG_MA + LCG_MB;
            sb = sb * LCG_MA + LCG_MB;
            as[i] = sa & mask;
            bs[i] = sb & mask;
        }
    }

    // REPEAT は wall time を REPEAT 倍するので TLE と相談して決める。
    constexpr int REPEAT = 1;
    uint64_t best_ns = ~uint64_t(0);
    u64 result_out = 0;

    for (int rep = 0; rep < REPEAT; ++rep) {
        u64 acc = 0;
        auto t0 = chrono::steady_clock::now();
        for (u64 i = 0; i < n_; ++i) {
            acc ^= G::gcd(as[i], bs[i]);
        }
        auto t1 = chrono::steady_clock::now();
        result_out = acc;
        auto ns = (uint64_t)chrono::duration_cast<chrono::nanoseconds>(t1 - t0).count();
        if (ns < best_ns) best_ns = ns;
    }

    fprintf(stderr, "ALGO_TIME_NS=%llu\n", (unsigned long long)best_ns);
    cout << result_out << '\n';
    return 0;
}
