#pragma once
// =============================================================================
// Source: yosupo "Pfaffian of Matrix" 提出 278726 (cp-algo).
//   https://judge.yosupo.jp/submission/278726
//   方式: skew-symmetric Gauss elimination で 2N×2N 交代行列を [[0,*][−*,0]]
//   ブロック対角に変形しながら、odd ステップごとに pivot 値を掛けて Pf を蓄積。
//   cp-algo は valarray-backed の vec<base> / matrix<base> で書かれており、
//   add_scaled で 4 回ごとに pseudonormalize の lazy reduce が入る。
//   この提出は古めの cp-algo 系 (FFT 未使用) なので依存連鎖が浅く、
//   modint / vector / matrix だけで完結する。
// 抽出方針:
//   - blazingio 等の I/O 高速化は使わず base.cpp 側の cin で済ませる
//   - #line マーカー除去
//   - Pfaffian::run wrapper で vector<vector<u32>> → matrix<base> に詰めて
//     元提出 solve() のロジックをそのまま実行
// ライセンス: cp-algo (元実装に従う)
// =============================================================================

#pragma GCC optimize("Ofast,unroll-loops")
#include "_common.hpp"







#include <chrono>
#include <random>
namespace cp_algo::random {
    uint64_t rng() {
        static std::mt19937_64 rng(
            std::chrono::steady_clock::now().time_since_epoch().count()
        );
        return rng();
    }
}




#include <functional>
#include <cstdint>
namespace cp_algo::math {
#ifdef CP_ALGO_MAXN
    const int maxn = CP_ALGO_MAXN;
#else
    const int maxn = 1 << 19;
#endif
    const int magic = 64; // threshold for sizes to run the naive algo

    auto bpow(auto const& x, auto n, auto const& one, auto op) {
        if(n == 0) {
            return one;
        } else {
            auto t = bpow(x, n / 2, one, op);
            t = op(t, t);
            if(n % 2) {
                t = op(t, x);
            }
            return t;
        }
    }
    auto bpow(auto x, auto n, auto ans) {
        return bpow(x, n, ans, std::multiplies{});
    }
    template<typename T>
    T bpow(T const& x, auto n) {
        return bpow(x, n, T(1));
    }
}








#include <iostream>
#include <cassert>
namespace cp_algo::math {
    inline constexpr auto inv2(auto x) {
        assert(x % 2);
        std::make_unsigned_t<decltype(x)> y = 1;
        while(y * x != 1) {
            y *= 2 - x * y;
        }
        return y;
    }

    template<typename modint, typename _Int>
    struct modint_base {
        using Int = _Int;
        using UInt = std::make_unsigned_t<Int>;
        static constexpr size_t bits = sizeof(Int) * 8;
        using Int2 = std::conditional_t<bits <= 32, int64_t, __int128_t>;
        using UInt2 = std::conditional_t<bits <= 32, uint64_t, __uint128_t>;
        static Int mod() {
            return modint::mod();
        }
        static UInt imod() {
            return modint::imod();
        }
        static UInt2 pw128() {
            return modint::pw128();
        }
        static UInt m_reduce(UInt2 ab) {
            if(mod() % 2 == 0) [[unlikely]] {
                return UInt(ab % mod());
            } else {
                UInt2 m = (UInt)ab * imod();
                return UInt((ab + m * mod()) >> bits);
            }
        }
        static UInt m_transform(UInt a) {
            if(mod() % 2 == 0) [[unlikely]] {
                return a;
            } else {
                return m_reduce(a * pw128());
            }
        }
        modint_base(): r(0) {}
        modint_base(Int2 rr): r(UInt(rr % mod())) {
            r = std::min(r, r + mod());
            r = m_transform(r);
        }
        modint inv() const {
            return bpow(to_modint(), mod() - 2);
        }
        modint operator - () const {
            modint neg;
            neg.r = std::min(-r, 2 * mod() - r);
            return neg;
        }
        modint& operator /= (const modint &t) {
            return to_modint() *= t.inv();
        }
        modint& operator *= (const modint &t) {
            r = m_reduce((UInt2)r * t.r);
            return to_modint();
        }
        modint& operator += (const modint &t) {
            r += t.r; r = std::min(r, r - 2 * mod());
            return to_modint();
        }
        modint& operator -= (const modint &t) {
            r -= t.r; r = std::min(r, r + 2 * mod());
            return to_modint();
        }
        modint operator + (const modint &t) const {return modint(to_modint()) += t;}
        modint operator - (const modint &t) const {return modint(to_modint()) -= t;}
        modint operator * (const modint &t) const {return modint(to_modint()) *= t;}
        modint operator / (const modint &t) const {return modint(to_modint()) /= t;}
        // Why <=> doesn't work?..
        auto operator == (const modint_base &t) const {return getr() == t.getr();}
        auto operator != (const modint_base &t) const {return getr() != t.getr();}
        auto operator <= (const modint_base &t) const {return getr() <= t.getr();}
        auto operator >= (const modint_base &t) const {return getr() >= t.getr();}
        auto operator < (const modint_base &t) const {return getr() < t.getr();}
        auto operator > (const modint_base &t) const {return getr() > t.getr();}
        Int rem() const {
            UInt R = getr();
            return 2 * R > (UInt)mod() ? R - mod() : R;
        }

        // Only use if you really know what you're doing!
        UInt modmod() const {return (UInt)8 * mod() * mod();};
        void add_unsafe(UInt t) {r += t;}
        void pseudonormalize() {r = std::min(r, r - modmod());}
        modint const& normalize() {
            if(r >= (UInt)mod()) {
                r %= mod();
            }
            return to_modint();
        }
        void setr(UInt rr) {r = m_transform(rr);}
        UInt getr() const {
            UInt res = m_reduce(r);
            return std::min(res, res - mod());
        }
        void setr_direct(UInt rr) {r = rr;}
        UInt getr_direct() const {return r;}
    private:
        UInt r;
        modint& to_modint() {return static_cast<modint&>(*this);}
        modint const& to_modint() const {return static_cast<modint const&>(*this);}
    };
    template<typename modint>
    concept modint_type = std::is_base_of_v<modint_base<modint, typename modint::Int>, modint>;
    template<modint_type modint>
    std::istream& operator >> (std::istream &in, modint &x) {
        typename modint::UInt r;
        auto &res = in >> r;
        x.setr(r);
        return res;
    }
    template<modint_type modint>
    std::ostream& operator << (std::ostream &out, modint const& x) {
        return out << x.getr();
    }

    template<auto m>
    struct modint: modint_base<modint<m>, decltype(m)> {
        using Base = modint_base<modint<m>, decltype(m)>;
        using Base::Base;
        static constexpr Base::UInt im = m % 2 ? inv2(-m) : 0;
        static constexpr Base::UInt r2 = (typename Base::UInt2)(-1) % m + 1;
        static constexpr Base::Int mod() {return m;}
        static constexpr Base::UInt imod() {return im;}
        static constexpr Base::UInt2 pw128() {return r2;}
    };

    template<typename Int = int64_t>
    struct dynamic_modint: modint_base<dynamic_modint<Int>, Int> {
        using Base = modint_base<dynamic_modint<Int>, Int>;
        using Base::Base;
        static Int mod() {return m;}
        static Base::UInt imod() {return im;}
        static Base::UInt2 pw128() {return r2;}
        static void switch_mod(Int nm) {
            m = nm;
            im = m % 2 ? inv2(-m) : 0;
            r2 = static_cast<Base::UInt>(static_cast<Base::UInt2>(-1) % m + 1);
        }

        // Wrapper for temp switching
        auto static with_mod(Int tmp, auto callback) {
            struct scoped {
                Int prev = mod();
                ~scoped() {switch_mod(prev);}
            } _;
            switch_mod(tmp);
            return callback();
        }
    private:
        static thread_local Int m;
        static thread_local Base::UInt im, r2;
    };
    template<typename Int>
    Int thread_local dynamic_modint<Int>::m = 1;
    template<typename Int>
    dynamic_modint<Int>::Base::UInt thread_local dynamic_modint<Int>::im = -1;
    template<typename Int>
    dynamic_modint<Int>::Base::UInt thread_local dynamic_modint<Int>::r2 = 0;
}


#include <algorithm>
#include <valarray>

#include <iterator>

namespace cp_algo::linalg {
    template<class vec, typename base>
    struct valarray_base: std::valarray<base> {
        using Base = std::valarray<base>;
        using Base::Base;

        valarray_base(base const& t): Base(t, 1) {}

        auto begin() {return std::begin(to_valarray());}
        auto begin() const {return std::begin(to_valarray());}
        auto end() {return std::end(to_valarray());}
        auto end() const {return std::end(to_valarray());}

        bool operator == (vec const& t) const {return std::ranges::equal(*this, t);}
        bool operator != (vec const& t) const {return !(*this == t);}

        vec operator-() const {return Base::operator-();}

        static vec from_range(auto const& R) {
            vec res(std::ranges::distance(R));
            std::ranges::copy(R, res.begin());
            return res;
        }
        Base& to_valarray() {return static_cast<Base&>(*this);}
        Base const& to_valarray() const {return static_cast<Base const&>(*this);}
    };

    template<class vec, typename base>
    vec operator+(valarray_base<vec, base> const& a, valarray_base<vec, base> const& b) {
        return a.to_valarray() + b.to_valarray();
    }
    template<class vec, typename base>
    vec operator-(valarray_base<vec, base> const& a, valarray_base<vec, base> const& b) {
        return a.to_valarray() - b.to_valarray();
    }

    template<class vec, typename base>
    struct vec_base: valarray_base<vec, base> {
        using Base = valarray_base<vec, base>;
        using Base::Base;

        static vec ei(size_t n, size_t i) {
            vec res(n);
            res[i] = 1;
            return res;
        }

        virtual void add_scaled(vec const& b, base scale, size_t i = 0) {
            if(scale != base(0)) {
                for(; i < size(*this); i++) {
                    (*this)[i] += scale * b[i];
                }
            }
        }
        virtual vec const& normalize() {
            return static_cast<vec&>(*this);
        }
        virtual base normalize(size_t i) {
            return (*this)[i];
        }
        void read() {
            for(auto &it: *this) {
                std::cin >> it;
            }
        }
        void print() const {
            std::ranges::copy(*this, std::ostream_iterator<base>(std::cout, " "));
            std::cout << "\n";
        }
        static vec random(size_t n) {
            vec res(n);
            std::ranges::generate(res, random::rng);
            return res;
        }
        // Concatenate vectors
        vec operator |(vec const& t) const {
            vec res(size(*this) + size(t));
            res[std::slice(0, size(*this), 1)] = *this;
            res[std::slice(size(*this), size(t), 1)] = t;
            return res;
        }

        // Generally, vec shouldn't be modified
        // after its pivot index is set
        std::pair<size_t, base> find_pivot() {
            if(pivot == size_t(-1)) {
                pivot = 0;
                while(pivot < size(*this) && normalize(pivot) == base(0)) {
                    pivot++;
                }
                if(pivot < size(*this)) {
                    pivot_inv = base(1) / (*this)[pivot];
                }
            }
            return {pivot, pivot_inv};
        }
        void reduce_by(vec &t) {
            auto [pivot, pinv] = t.find_pivot();
            if(pivot < size(*this)) {
                add_scaled(t, -normalize(pivot) * pinv, pivot);
            }
        }
    private:
        size_t pivot = -1;
        base pivot_inv;
    };

    template<typename base>
    struct vec: vec_base<vec<base>, base> {
        using Base = vec_base<vec<base>, base>;
        using Base::Base;
    };

    template<math::modint_type base>
    struct vec<base>: vec_base<vec<base>, base> {
        using Base = vec_base<vec<base>, base>;
        using Base::Base;

        void add_scaled(vec const& b, base scale, size_t i = 0) override {
            static_assert(base::bits >= 64, "Only wide modint types for linalg");
            uint64_t scaler = scale.getr();
            if(scale != base(0)) {
                for(; i < size(*this); i++) {
                    (*this)[i].add_unsafe(scaler * b[i].getr_direct());
                }
                if(++counter == 4) {
                    for(auto &it: *this) {
                        it.pseudonormalize();
                    }
                    counter = 0;
                }
            }
        }
        vec const& normalize() override {
            for(auto &it: *this) {
                it.normalize();
            }
            return *this;
        }
        base normalize(size_t i) override {
            return (*this)[i].normalize();
        }
    private:
        size_t counter = 0;
    };
}


#include <optional>

#include <vector>
#include <array>
namespace cp_algo::linalg {
    enum gauss_mode {normal, reverse};
    template<typename base_t>
    struct matrix: valarray_base<matrix<base_t>, vec<base_t>> {
        using base = base_t;
        using Base = valarray_base<matrix<base>, vec<base>>;
        using Base::Base;

        matrix(size_t n): Base(vec<base>(n), n) {}
        matrix(size_t n, size_t m): Base(vec<base>(m), n) {}

        size_t n() const {return size(*this);}
        size_t m() const {return n() ? size(row(0)) : 0;}
        auto dim() const {return std::array{n(), m()};}

        auto& row(size_t i) {return (*this)[i];}
        auto const& row(size_t i) const {return (*this)[i];}

        matrix& operator *=(base t) {for(auto &it: *this) it *= t; return *this;}
        matrix operator *(base t) const {return matrix(*this) *= t;}
        matrix& operator /=(base t) {return *this *= base(1) / t;}
        matrix operator /(base t) const {return matrix(*this) /= t;}

        // Make sure the result is matrix, not Base
        matrix& operator *=(matrix const& t) {return *this = *this * t;}

        void read() {
            for(auto &it: *this) {
                it.read();
            }
        }
        void print() const {
            for(auto const& it: *this) {
                it.print();
            }
        }

        static matrix block_diagonal(std::vector<matrix> const& blocks) {
            size_t n = 0;
            for(auto &it: blocks) {
                assert(it.n() == it.m());
                n += it.n();
            }
            matrix res(n);
            n = 0;
            for(auto &it: blocks) {
                for(size_t i = 0; i < it.n(); i++) {
                    res[n + i][std::slice(n, it.n(), 1)] = it[i];
                }
                n += it.n();
            }
            return res;
        }
        static matrix random(size_t n, size_t m) {
            matrix res(n, m);
            std::ranges::generate(res, std::bind(vec<base>::random, m));
            return res;
        }
        static matrix random(size_t n) {
            return random(n, n);
        }
        static matrix eye(size_t n) {
            matrix res(n);
            for(size_t i = 0; i < n; i++) {
                res[i][i] = 1;
            }
            return res;
        }

        // Concatenate matrices
        matrix operator |(matrix const& b) const {
            assert(n() == b.n());
            matrix res(n(), m()+b.m());
            for(size_t i = 0; i < n(); i++) {
                res[i] = row(i) | b[i];
            }
            return res;
        }
        matrix submatrix(auto slicex, auto slicey) const {
            matrix res = (*this)[slicex];
            for(auto &row: res) {
                row = vec<base>(row[slicey]);
            }
            return res;
        }

        matrix T() const {
            matrix res(m(), n());
            for(size_t i = 0; i < n(); i++) {
                for(size_t j = 0; j < m(); j++) {
                    res[j][i] = row(i)[j];
                }
            }
            return res;
        }

        matrix operator *(matrix const& b) const {
            assert(m() == b.n());
            matrix res(n(), b.m());
            for(size_t i = 0; i < n(); i++) {
                for(size_t j = 0; j < m(); j++) {
                    res[i].add_scaled(b[j], row(i)[j]);
                }
            }
            return res.normalize();
        }

        vec<base> apply(vec<base> const& x) const {
            return (matrix(x) * *this)[0];
        }

        matrix pow(uint64_t k) const {
            assert(n() == m());
            return bpow(*this, k, eye(n()));
        }

        matrix& normalize() {
            for(auto &it: *this) {
                it.normalize();
            }
            return *this;
        }
        template<gauss_mode mode = normal>
        void eliminate(size_t i, size_t k) {
            auto kinv = base(1) / row(i).normalize()[k];
            for(size_t j = (mode == normal) * i; j < n(); j++) {
                if(j != i) {
                    row(j).add_scaled(row(i), -row(j).normalize(k) * kinv);
                }
            }
        }
        template<gauss_mode mode = normal>
        void eliminate(size_t i) {
            row(i).normalize();
            for(size_t j = (mode == normal) * i; j < n(); j++) {
                if(j != i) {
                    row(j).reduce_by(row(i));
                }
            }
        }
        template<gauss_mode mode = normal>
        matrix& gauss() {
            for(size_t i = 0; i < n(); i++) {
                eliminate<mode>(i);
            }
            return normalize();
        }
        template<gauss_mode mode = normal>
        auto echelonize(size_t lim) {
            return gauss<mode>().sort_classify(lim);
        }
        template<gauss_mode mode = normal>
        auto echelonize() {
            return echelonize<mode>(m());
        }

        size_t rank() const {
            if(n() > m()) {
                return T().rank();
            }
            return size(matrix(*this).echelonize()[0]);
        }

        base det() const {
            assert(n() == m());
            matrix b = *this;
            b.echelonize();
            base res = 1;
            for(size_t i = 0; i < n(); i++) {
                res *= b[i][i];
            }
            return res;
        }

        std::pair<base, matrix> inv() const {
            assert(n() == m());
            matrix b = *this | eye(n());
            if(size(b.echelonize<reverse>(n())[0]) < n()) {
                return {0, {}};
            }
            base det = 1;
            for(size_t i = 0; i < n(); i++) {
                det *= b[i][i];
                b[i] *= base(1) / b[i][i];
            }
            return {det, b.submatrix(std::slice(0, n(), 1), std::slice(n(), n(), 1))};
        }

        // Can also just run gauss on T() | eye(m)
        // but it would be slower :(
        auto kernel() const {
            auto A = *this;
            auto [pivots, free] = A.template echelonize<reverse>();
            matrix sols(size(free), m());
            for(size_t j = 0; j < size(pivots); j++) {
                base scale = base(1) / A[j][pivots[j]];
                for(size_t i = 0; i < size(free); i++) {
                    sols[i][pivots[j]] = A[j][free[i]] * scale;
                }
            }
            for(size_t i = 0; i < size(free); i++) {
                sols[i][free[i]] = -1;
            }
            return sols;
        }

        // [solution, basis], transposed
        std::optional<std::array<matrix, 2>> solve(matrix t) const {
            matrix sols = (*this | t).kernel();
            if(sols.n() < t.m() || sols.submatrix(
                std::slice(sols.n() - t.m(), t.m(), 1),
                std::slice(m(), t.m(), 1)
            ) != -eye(t.m())) {
                return std::nullopt;
            } else {
                return std::array{
                    sols.submatrix(std::slice(sols.n() - t.m(), t.m(), 1),
                                   std::slice(0, m(), 1)),
                    sols.submatrix(std::slice(0, sols.n() - t.m(), 1),
                                   std::slice(0, m(), 1))
                };
            }
        }

        // To be called after a gaussian elimination run
        // Sorts rows by pivots and classifies
        // variables into pivots and free
        auto sort_classify(size_t lim) {
            size_t rk = 0;
            std::vector<size_t> free, pivots;
            for(size_t j = 0; j < lim; j++) {
                for(size_t i = rk + 1; i < n() && row(rk)[j] == base(0); i++) {
                    if(row(i)[j] != base(0)) {
                        std::swap(row(i), row(rk));
                        row(rk) = -row(rk);
                    }
                }
                if(rk < n() && row(rk)[j] != base(0)) {
                    pivots.push_back(j);
                    rk++;
                } else {
                    free.push_back(j);
                }
            }
            return std::array{pivots, free};
        }
    };
    template<typename base_t>
    auto operator *(base_t t, matrix<base_t> const& A) {return A * t;}
}


#include <bits/stdc++.h>

using namespace std;
using namespace cp_algo;

const int64_t mod = 998244353;
using base = math::modint<mod>;

// ---------------------------------------------------------------------------
// Wrapper
// ---------------------------------------------------------------------------
struct Pfaffian {
 static u32 run(int N, const vector<vector<u32>>& M_in) {
  using base = cp_algo::math::modint<(int64_t) MOD>;
  int n = 2 * N;
  cp_algo::linalg::matrix<base> A((size_t) n);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j) A[i][j].setr((uint64_t) M_in[i][j]);
  base res = 1;
  for (int i = 1; i < n; ++i) {
   for (int j = i + 1; j < n; ++j) {
    if (A[j].normalize(i - 1) != base(0)) {
     std::swap(A[i], A[j]);
     for (int k = i; k < n; ++k) std::swap(A[k][i], A[k][j]);
     res *= base(-1);
     break;
    }
   }
   A[i].normalize();
   if (i % 2 == 1) res *= -A[i][i - 1];
   if (A[i][i - 1] == base(0)) return 0;
   base Ai = A[i][i - 1].inv();
   for (int j = i + 1; j < n; ++j) {
    A[j].add_scaled(A[i], -A[j].normalize(i - 1) * Ai, i);
   }
  }
  u32 r = (u32) res.getr();
  if (r >= MOD) r -= MOD;
  return r;
 }
};
