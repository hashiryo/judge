// harness: 各 algos/*.hpp が定義する struct MP を使って pow(a_i, b_i) を計測する。
// runtime-62 (64-bit mod) の 30-bit 素数版。a_i, b_i ともに u32 で持つ。
#include "algos/_common.hpp"

#ifndef ALGO_HPP
#define ALGO_HPP "algos/naive32.hpp"
#endif
#include ALGO_HPP

signed main() {
    cin.tie(0);
    ios::sync_with_stdio(false);
    u64 n_, mod_, bbits_, sa_, sb_;
    cin >> n_ >> mod_ >> bbits_ >> sa_ >> sb_;

    constexpr u64 LCG_MA = 6364136223846793005ULL;
    constexpr u64 LCG_MB = 1442695040888963407ULL;
    const u32 b_mask = bbits_ >= 32 ? ~u32(0) : ((u32(1) << bbits_) - 1);
    vector<u32> as(n_), bs(n_);
    {
        u64 sa = sa_, sb = sb_;
        for (u64 i = 0; i < n_; ++i) {
            sa = sa * LCG_MA + LCG_MB;
            sb = sb * LCG_MA + LCG_MB;
            as[i] = u32(sa % mod_);
            bs[i] = u32(sb) & b_mask;
        }
    }

    constexpr int REPEAT = 1;
    uint64_t best_ns = ~uint64_t(0);
    u32 result_out = 0;

    for (int rep = 0; rep < REPEAT; ++rep) {
        u32 acc = 0;
        auto t0 = chrono::steady_clock::now();
        const MP mp{u32(mod_)};  // precompute も計測対象
        for (u64 i = 0; i < n_; ++i) {
            acc ^= mp.get(mp.pow(mp.set(as[i]), bs[i]));
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
