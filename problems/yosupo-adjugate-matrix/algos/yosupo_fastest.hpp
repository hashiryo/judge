#pragma once
// =============================================================================
// Source: yosupo "Adjugate Matrix" 提出 285779 (cp-algo).
//   https://judge.yosupo.jp/submission/285779
//   方式: A を (n+1)×(n+1) に乱数で拡張し、その逆行列 Ai と det D を取って
//   adj(A)[i][j] = D * (Ai[n][n] * Ai[i][j] - Ai[i][n] * Ai[n][j])。
//   D == 0 (拡張行列も特異) の場合は全 0 を出すが、確率的にほぼ起きない。
//   matrix.hpp は system 版 (#285778) と同じ。modint_vec の add_scaled で
//   AVX2 (_mm256_mul_epu32) montgomery 風 mul-add、4 回ごとに lazy reduce。
// 抽出方針:
//   - blazingio I/O は使わない (base.cpp が cin で済ませる)
//   - <experimental/simd> include を削除 (libstdc++ stdfloat::float16_t と
//     arm_neon.h ::float16_t の衝突回避; cp-algo 自身は stdx::* を使っていない)
//   - #line マーカー除去
//   - Rank::run wrapper で vector<vector<u32>> → matrix<base> に詰める
//     (base::setr で u32 をそのまま入れる)
// ライセンス: cp-algo (元実装に従う)
// =============================================================================

#pragma GCC optimize("Ofast,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2,bmi,bmi2")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#elif defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "_common.hpp"

#include <ranges>








namespace cp_algo::random {
    uint64_t rng() {
        static std::mt19937_64 rng(
            std::chrono::steady_clock::now().time_since_epoch().count()
        );
        return rng();
    }
}





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








namespace cp_algo::math {

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
        static Int remod() {
            return modint::remod();
        }
        static UInt2 modmod() {
            return UInt2(mod()) * mod();
        }
        modint_base() = default;
        modint_base(Int2 rr) {
            to_modint().setr(UInt((rr + modmod()) % mod()));
        }
        modint inv() const {
            return bpow(to_modint(), mod() - 2);
        }
        modint operator - () const {
            modint neg;
            neg.r = std::min(-r, remod() - r);
            return neg;
        }
        modint& operator /= (const modint &t) {
            return to_modint() *= t.inv();
        }
        modint& operator *= (const modint &t) {
            r = UInt(UInt2(r) * t.r % mod());
            return to_modint();
        }
        modint& operator += (const modint &t) {
            r += t.r; r = std::min(r, r - remod());
            return to_modint();
        }
        modint& operator -= (const modint &t) {
            r -= t.r; r = std::min(r, r + remod());
            return to_modint();
        }
        modint operator + (const modint &t) const {return modint(to_modint()) += t;}
        modint operator - (const modint &t) const {return modint(to_modint()) -= t;}
        modint operator * (const modint &t) const {return modint(to_modint()) *= t;}
        modint operator / (const modint &t) const {return modint(to_modint()) /= t;}
        // Why <=> doesn't work?..
        auto operator == (const modint &t) const {return to_modint().getr() == t.getr();}
        auto operator != (const modint &t) const {return to_modint().getr() != t.getr();}
        auto operator <= (const modint &t) const {return to_modint().getr() <= t.getr();}
        auto operator >= (const modint &t) const {return to_modint().getr() >= t.getr();}
        auto operator < (const modint &t) const {return to_modint().getr() < t.getr();}
        auto operator > (const modint &t) const {return to_modint().getr() > t.getr();}
        Int rem() const {
            UInt R = to_modint().getr();
            return R - (R > (UInt)mod() / 2) * mod();
        }
        void setr(UInt rr) {
            r = rr;
        }
        UInt getr() const {
            return r;
        }

        // Only use these if you really know what you're doing!
        static UInt modmod8() {return UInt(8 * modmod());}
        void add_unsafe(UInt t) {r += t;}
        void pseudonormalize() {r = std::min(r, r - modmod8());}
        modint const& normalize() {
            if(r >= (UInt)mod()) {
                r %= mod();
            }
            return to_modint();
        }
        void setr_direct(UInt rr) {r = rr;}
        UInt getr_direct() const {return r;}
    protected:
        UInt r;
    private:
        modint& to_modint() {return static_cast<modint&>(*this);}
        modint const& to_modint() const {return static_cast<modint const&>(*this);}
    };
    template<typename modint>
    concept modint_type = std::is_base_of_v<modint_base<modint, typename modint::Int>, modint>;
    template<modint_type modint>
    decltype(std::cin)& operator >> (decltype(std::cin) &in, modint &x) {
        typename modint::UInt r;
        auto &res = in >> r;
        x.setr(r);
        return res;
    }
    template<modint_type modint>
    decltype(std::cout)& operator << (decltype(std::cout) &out, modint const& x) {
        return out << x.getr();
    }

    template<auto m>
    struct modint: modint_base<modint<m>, decltype(m)> {
        using Base = modint_base<modint<m>, decltype(m)>;
        using Base::Base;
        static constexpr Base::Int mod() {return m;}
        static constexpr Base::UInt remod() {return m;}
        auto getr() const {return Base::r;}
    };

    inline constexpr auto inv2(auto x) {
        assert(x % 2);
        std::make_unsigned_t<decltype(x)> y = 1;
        while(y * x != 1) {
            y *= 2 - x * y;
        }
        return y;
    }

    template<typename Int = int64_t>
    struct dynamic_modint: modint_base<dynamic_modint<Int>, Int> {
        using Base = modint_base<dynamic_modint<Int>, Int>;
        using Base::Base;

        static Base::UInt m_reduce(Base::UInt2 ab) {
            if(mod() % 2 == 0) [[unlikely]] {
                return typename Base::UInt(ab % mod());
            } else {
                typename Base::UInt2 m = typename Base::UInt(ab) * imod();
                return typename Base::UInt((ab + m * mod()) >> Base::bits);
            }
        }
        static Base::UInt m_transform(Base::UInt a) {
            if(mod() % 2 == 0) [[unlikely]] {
                return a;
            } else {
                return m_reduce(a * pw128());
            }
        }
        dynamic_modint& operator *= (const dynamic_modint &t) {
            Base::r = m_reduce(typename Base::UInt2(Base::r) * t.r);
            return *this;
        }
        void setr(Base::UInt rr) {
            Base::r = m_transform(rr);
        }
        Base::UInt getr() const {
            typename Base::UInt res = m_reduce(Base::r);
            return std::min(res, res - mod());
        }
        static Int mod() {return m;}
        static Int remod() {return 2 * m;}
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







// Single macro to detect POSIX platforms (Linux, Unix, macOS)
// macOS は MADV_HUGEPAGE / MADV_POPULATE_WRITE がないので mmap path 無効化
#if defined(__linux__)
#  define CP_ALGO_USE_MMAP 1
#  include <sys/mman.h>
#else
#  define CP_ALGO_USE_MMAP 0
#endif

namespace cp_algo {
    template <typename T, std::size_t Align = 32>
    class big_alloc {
        static_assert( Align >= alignof(void*), "Align must be at least pointer-size");
        static_assert(std::popcount(Align) == 1, "Align must be a power of two");
    public:
        using value_type = T;
        template <class U> struct rebind { using other = big_alloc<U, Align>; };
        constexpr bool operator==(const big_alloc&) const = default;
        constexpr bool operator!=(const big_alloc&) const = default;

        big_alloc() noexcept = default;
        template <typename U, std::size_t A>
        big_alloc(const big_alloc<U, A>&) noexcept {}

        [[nodiscard]] T* allocate(std::size_t n) {
            std::size_t padded = round_up(n * sizeof(T));
            std::size_t align = std::max<std::size_t>(alignof(T),  Align);
#if CP_ALGO_USE_MMAP
            if (padded >= MEGABYTE) {
                void* raw = mmap(nullptr, padded,
                                PROT_READ | PROT_WRITE,
                                MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
                madvise(raw, padded, MADV_HUGEPAGE);
                madvise(raw, padded, MADV_POPULATE_WRITE);
                return static_cast<T*>(raw);
            }
#endif
            return static_cast<T*>(::operator new(padded, std::align_val_t(align)));
        }

        void deallocate(T* p, std::size_t n) noexcept {
            if (!p) return;
            std::size_t padded = round_up(n * sizeof(T));
            std::size_t align  = std::max<std::size_t>(alignof(T),  Align);
    #if CP_ALGO_USE_MMAP
            if (padded >= MEGABYTE) { munmap(p, padded); return; }
    #endif
            ::operator delete(p, padded, std::align_val_t(align));
        }

    private:
        static constexpr std::size_t MEGABYTE = 1 << 20;
        static constexpr std::size_t round_up(std::size_t x) noexcept {
            return (x + Align - 1) / Align * Align;
        }
    };
}






namespace cp_algo {
    template<typename T, size_t len>
    using simd [[gnu::vector_size(len * sizeof(T))]] = T;
    using i64x4 = simd<int64_t, 4>;
    using u64x4 = simd<uint64_t, 4>;
    using u32x8 = simd<uint32_t, 8>;
    using i32x4 = simd<int32_t, 4>;
    using u32x4 = simd<uint32_t, 4>;
    using dx4 = simd<double, 4>;

    [[gnu::always_inline]] inline dx4 abs(dx4 a) {
    return a < 0 ? -a : a;
    }

    // https://stackoverflow.com/a/77376595
    // works for ints in (-2^51, 2^51)
    static constexpr dx4 magic = dx4() + (3ULL << 51);
    [[gnu::always_inline]] inline i64x4 lround(dx4 x) {
        return i64x4(x + magic) - i64x4(magic);
    }
    [[gnu::always_inline]] inline dx4 to_double(i64x4 x) {
        return dx4(x + i64x4(magic)) - magic;
    }

    [[gnu::always_inline]] inline dx4 round(dx4 a) {
        return dx4{
            std::nearbyint(a[0]),
            std::nearbyint(a[1]),
            std::nearbyint(a[2]),
            std::nearbyint(a[3])
        };
    }

    [[gnu::always_inline]] inline u64x4 montgomery_reduce(u64x4 x, u64x4 mod, u64x4 imod) {
        auto x_ninv = u64x4(u32x8(x) * u32x8(imod));
#ifdef __AVX2__
        x += u64x4(_mm256_mul_epu32(__m256i(x_ninv), __m256i(mod)));
#else
        x += x_ninv * mod;
#endif
        return x >> 32;
    }

    [[gnu::always_inline]] inline u64x4 montgomery_mul(u64x4 x, u64x4 y, u64x4 mod, u64x4 imod) {
#ifdef __AVX2__
        return montgomery_reduce(u64x4(_mm256_mul_epu32(__m256i(x), __m256i(y))), mod, imod);
#else
        return montgomery_reduce(x * y, mod, imod);
#endif
    }

    [[gnu::always_inline]] inline dx4 rotate_right(dx4 x) {
#if defined(__clang__)
        return __builtin_shufflevector(x, x, 3, 0, 1, 2);
#else
        static constexpr u64x4 shuffler = {3, 0, 1, 2};
        return __builtin_shuffle(x, shuffler);
#endif
    }

    template<class Target>
    [[gnu::always_inline]] inline Target& vector_cast(auto &&p) {
        return *reinterpret_cast<Target*>(std::assume_aligned<alignof(Target)>(&p));
    }
}





namespace cp_algo {
    std::map<std::string, double> checkpoints;
    template<bool final = false>
    void checkpoint([[maybe_unused]] std::string const& msg = "") {
#ifdef CP_ALGO_CHECKPOINT
        static double last = 0;
        double now = (double)clock() / CLOCKS_PER_SEC;
        double delta = now - last;
        last = now;
        if(msg.size() && !final) {
            checkpoints[msg] += delta;
        }
        if(final) {
            for(auto const& [key, value] : checkpoints) {
                std::cerr << key << ": " << value * 1000 << " ms\n";
            }
            std::cerr << "Total: " << now * 1000 << " ms\n";
        }
#endif
    }
}


#include <ranges>
namespace cp_algo::linalg {
    template<typename base, class Alloc = big_alloc<base>>
    struct vec: std::basic_string<base, std::char_traits<base>, Alloc> {
        using Base = std::basic_string<base, std::char_traits<base>, Alloc>;
        using Base::Base;

        vec(Base const& t): Base(t) {}
        vec(Base &&t): Base(std::move(t)) {}
        vec(size_t n): Base(n, base()) {}
        vec(auto &&r): Base(std::ranges::to<Base>(r)) {}

        static vec ei(size_t n, size_t i) {
            vec res(n);
            res[i] = 1;
            return res;
        }

        auto operator-() const {
            return *this | std::views::transform([](auto x) {return -x;});
        }
        auto operator *(base t) const {
            return *this | std::views::transform([t](auto x) {return x * t;});
        }
        auto operator *=(base t) {
            for(auto &it: *this) {
                it *= t;
            }
            return *this;
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
            for(auto &it: *this) {
                std::cout << it << " ";
            }
            std::cout << "\n";
        }
        static vec random(size_t n) {
            vec res(n);
            std::ranges::generate(res, random::rng);
            return res;
        }
        // Concatenate vectors
        vec operator |(vec const& t) const {
            return std::views::join(std::array{
                std::views::all(*this),
                std::views::all(t)
            });
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

    template<math::modint_type base, class Alloc = big_alloc<base>>
    struct modint_vec: vec<base, Alloc> {
        using Base = vec<base, Alloc>;
        using Base::Base;

        modint_vec(Base const& t): Base(t) {}
        modint_vec(Base &&t): Base(std::move(t)) {}

        void add_scaled(Base const& b, base scale, size_t i = 0) override {
            static_assert(base::bits >= 64, "Only wide modint types for linalg");
            if(scale != base(0)) {
                assert(Base::size() == b.size());
                size_t n = size(*this);
                u64x4 scaler = u64x4() + scale.getr();
                for(i -= i % 4; i < n; i += 4) {
                    auto &ai = vector_cast<u64x4>((*this)[i]);
                    auto bi = vector_cast<u64x4 const>(b[i]);
#ifdef __AVX2__
                    ai += u64x4(_mm256_mul_epu32(__m256i(scaler), __m256i(bi)));
#else
                    ai += scaler * bi;
#endif
                }
                if(++counter == 4) {
                    for(auto &it: *this) {
                        it.pseudonormalize();
                    }
                    counter = 0;
                }
            }
        }
        Base const& normalize() override {
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

namespace cp_algo::linalg {
    enum gauss_mode {normal, reverse};

    template<typename base_t, class _vec_t = std::conditional_t<
        math::modint_type<base_t>,
        modint_vec<base_t>,
        vec<base_t>>>
    struct matrix: std::vector<_vec_t> {
        using vec_t = _vec_t;
        using base = base_t;
        using Base = std::vector<vec_t>;
        using Base::Base;

        matrix(size_t n): Base(n, vec_t(n)) {}
        matrix(size_t n, size_t m): Base(n, vec_t(m)) {}

        matrix(Base const& t): Base(t) {}
        matrix(Base &&t): Base(std::move(t)) {}

        static matrix from(auto &&r) {
            return std::ranges::to<Base>(r);
        }

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

        void read_transposed() {
            for(size_t j = 0; j < m(); j++) {
                for(size_t i = 0; i < n(); i++) {
                    std::cin >> (*this)[i][j];
                }
            }
        }
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
                    std::ranges::copy(it[i], begin(res[n + i]) + n);
                }
                n += it.n();
            }
            return res;
        }
        static matrix random(size_t n, size_t m) {
            matrix res(n, m);
            std::ranges::generate(res, std::bind(vec_t::random, m));
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
        // 元提出 (system) の submatrix: rows と cols の views を 1 段で適用。
        matrix submatrix(auto viewx, auto viewy) const {
            return from(*this | viewx | std::views::transform(
                [&](auto const& y) {return vec_t(y | viewy);}));
        }
        // 元提出 (system) には matrix の単項マイナスがある。-eye(...) で使う。
        auto operator-() const {
            return from(*this | std::views::transform([](auto x) {return vec_t(-x);}));
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

        vec_t apply(vec_t const& x) const {
            return (matrix(1, x) * *this)[0];
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
            return {det, b.submatrix(std::views::take(n()), std::views::drop(n()) | std::views::take(n()))};
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

        // [solution, basis], transposed (元提出 system の版)
        std::optional<std::array<matrix, 2>> solve(matrix t) const {
            matrix sols = (*this | t).kernel();
            if(sols.n() < t.m() || sols.submatrix(
                std::views::drop(sols.n() - t.m()),
                std::views::drop(m())
            ) != -eye(t.m())) {
                return std::nullopt;
            } else {
                return std::array{
                    sols.submatrix(std::views::drop(sols.n() - t.m()),
                                   std::views::take(m())),
                    sols.submatrix(std::views::take(sols.n() - t.m()),
                                   std::views::take(m()))
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

// ---------------------------------------------------------------------------
// Wrapper
// ---------------------------------------------------------------------------
struct Adjugate {
 static vector<vector<u32>> run(int n, const vector<vector<u32>>& a_in) {
  using base = cp_algo::math::modint<(int64_t) MOD>;
  // (n+1)×(n+1) に乱数で拡張
  cp_algo::linalg::matrix<base> A(n + 1, n + 1);
  for (int i = 0; i < n; ++i)
   for (int j = 0; j < n; ++j) A[i][j].setr((uint64_t) a_in[i][j]);
  for (int i = 0; i < n; ++i) {
   A[i][n].setr(cp_algo::random::rng() % (uint64_t) MOD);
   A[n][i].setr(cp_algo::random::rng() % (uint64_t) MOD);
  }
  A[n][n].setr(cp_algo::random::rng() % (uint64_t) MOD);
  auto [D, Ai] = A.inv();
  vector<vector<u32>> ans(n, vector<u32>(n, 0));
  if (D == base(0)) return ans; // ほぼ起きない
  for (int i = 0; i < n; ++i) {
   for (int j = 0; j < n; ++j) {
    base v = (Ai[n][n] * Ai[i][j] - Ai[i][n] * Ai[n][j]) * D;
    u32 r = (u32) v.getr();
    if (r >= MOD) r -= MOD;
    ans[i][j] = r;
   }
  }
  return ans;
 }
};
