#pragma once
// =============================================================================
// Source: yosupo "Exp of Set Power Series" 提出 339836 (cp-algo).
//   https://judge.yosupo.jp/submission/339836
//   方式: cp-algo の math::subset_exp<base>。subset_convolution と同じ
//   AVX2 montgomery rank-vector kernel を半分ずつの再帰で呼ぶ
//     out0 = subset_exp(g[:N/2])
//     out1 = subset_convolution(out0, g[N/2:])
//     return concat(out0, out1)
//   matrix.hpp 部分は subset_convolution 版 (#352953) より新しく、
//   xor_transform に transform_dir (forw/inv) enum が入っている。
// 抽出方針 (subset_convolution と同じ):
//   - blazingio / FastIO __io は使わず base.cpp 側の cin で済ませる
//   - <experimental/simd> 削除、mmap path を Linux 限定、__builtin_shuffle 分岐
//   - <generator> を __has_include で gate
//   - 未使用の subset_compose / subset_power_projection / read_bits 系を削除
//   - kth_set_bit (BMI2 _pdep_u64) を __BMI2__ で gate
// ライセンス: cp-algo (元実装に従う)
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#if (defined(__x86_64__) || defined(__i386__)) && !defined(USE_SIMDE)
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#endif
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#elif defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "_common.hpp"

#include <span>
#include <ranges>







#include <functional>
#include <cstdint>
#include <cassert>
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
    inline constexpr auto inv2(auto x) {
        assert(x % 2);
        std::make_unsigned_t<decltype(x)> y = 1;
        while(y * x != 1) {
            y *= 2 - x * y;
        }
        return y;
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
        constexpr static Int mod() {
            return modint::mod();
        }
        constexpr static Int remod() {
            return modint::remod();
        }
        constexpr static UInt2 modmod() {
            return UInt2(mod()) * mod();
        }
        constexpr modint_base() = default;
        constexpr modint_base(Int2 rr) {
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
        constexpr void setr(UInt rr) {
            r = rr;
        }
        constexpr UInt getr() const {
            return r;
        }

        // Only use these if you really know what you're doing!
        static uint64_t modmod8() {return uint64_t(8 * modmod());}
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
        constexpr modint& to_modint() {return static_cast<modint&>(*this);}
        constexpr modint const& to_modint() const {return static_cast<modint const&>(*this);}
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

    template<typename Int = int>
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









#include <cstddef>
#include <memory>

#if defined(__x86_64__) && !defined(CP_ALGO_DISABLE_AVX2) && !defined(USE_SIMDE)
#define CP_ALGO_SIMD_AVX2_TARGET _Pragma("GCC target(\"avx2\")")
#else
#define CP_ALGO_SIMD_AVX2_TARGET
#endif

#define CP_ALGO_SIMD_PRAGMA_PUSH \
    _Pragma("GCC push_options") \
    CP_ALGO_SIMD_AVX2_TARGET

CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo {
    template<typename T, size_t len>
    using simd [[gnu::vector_size(len * sizeof(T))]] = T;
    using u64x8 = simd<uint64_t, 8>;
    using u32x16 = simd<uint32_t, 16>;
    using i64x4 = simd<int64_t, 4>;
    using u64x4 = simd<uint64_t, 4>;
    using u32x8 = simd<uint32_t, 8>;
    using u16x16 = simd<uint16_t, 16>;
    using i32x4 = simd<int32_t, 4>;
    using u32x4 = simd<uint32_t, 4>;
    using u16x8 = simd<uint16_t, 8>;
    using u16x4 = simd<uint16_t, 4>;
    using i16x4 = simd<int16_t, 4>;
    using u8x32 = simd<uint8_t, 32>;
    using u8x8 = simd<uint8_t, 8>;
    using u8x4 = simd<uint8_t, 4>;
    using dx4 = simd<double, 4>;

    inline dx4 abs(dx4 a) {
        return dx4{
            std::abs(a[0]),
            std::abs(a[1]),
            std::abs(a[2]),
            std::abs(a[3])
        };
    }

    // https://stackoverflow.com/a/77376595
    // works for ints in (-2^51, 2^51)
    static constexpr dx4 magic = dx4() + (3ULL << 51);
    inline i64x4 lround(dx4 x) {
        return i64x4(x + magic) - i64x4(magic);
    }
    inline dx4 to_double(i64x4 x) {
        return dx4(x + i64x4(magic)) - magic;
    }

    inline dx4 round(dx4 a) {
        return dx4{
            std::nearbyint(a[0]),
            std::nearbyint(a[1]),
            std::nearbyint(a[2]),
            std::nearbyint(a[3])
        };
    }

    inline u64x4 low32(u64x4 x) {
        return x & uint32_t(-1);
    }
    inline auto swap_bytes(auto x) {
        return decltype(x)(__builtin_shufflevector(u32x8(x), u32x8(x), 1, 0, 3, 2, 5, 4, 7, 6));
    }
    inline u64x4 montgomery_reduce(u64x4 x, uint32_t mod, uint32_t imod) {
#ifdef __AVX2__
        auto x_ninv = u64x4(_mm256_mul_epu32(__m256i(x), __m256i() + imod));
        x += u64x4(_mm256_mul_epu32(__m256i(x_ninv), __m256i() + mod));
#else
        auto x_ninv = u64x4(u32x8(low32(x)) * imod);
        x += x_ninv * uint64_t(mod);
#endif
        return swap_bytes(x);
    }

    inline u64x4 montgomery_mul(u64x4 x, u64x4 y, uint32_t mod, uint32_t imod) {
#ifdef __AVX2__
        return montgomery_reduce(u64x4(_mm256_mul_epu32(__m256i(x), __m256i(y))), mod, imod);
#else
        return montgomery_reduce(x * y, mod, imod);
#endif
    }
    inline u32x8 montgomery_mul(u32x8 x, u32x8 y, uint32_t mod, uint32_t imod) {
        return u32x8(montgomery_mul(u64x4(x), u64x4(y), mod, imod)) |
               u32x8(swap_bytes(montgomery_mul(u64x4(swap_bytes(x)), u64x4(swap_bytes(y)), mod, imod)));
    }
    inline dx4 rotate_right(dx4 x) {
#if defined(__clang__)
        return __builtin_shufflevector(x, x, 3, 0, 1, 2);
#else
        static constexpr u64x4 shuffler = {3, 0, 1, 2};
        return __builtin_shuffle(x, shuffler);
#endif
    }

    template<std::size_t Align = 32>
    inline bool is_aligned(const auto* p) noexcept {
        return (reinterpret_cast<std::uintptr_t>(p) % Align) == 0;
    }

    template<class Target>
    inline Target& vector_cast(auto &&p) {
        return *reinterpret_cast<Target*>(std::assume_aligned<alignof(Target)>(&p));
    }
}
#pragma GCC pop_options





#include <set>
#include <map>
#include <deque>
#include <stack>
#include <queue>
#include <vector>
#include <string>

// libc++ (macOS) には <generator> がまだ無いので __has_include で gate。
#if __has_include(<generator>)
#include <generator>
#define CP_ALGO_HAS_GENERATOR 1
#endif
#include <forward_list>

// macOS は MADV_HUGEPAGE / MADV_POPULATE_WRITE が無いので Linux 限定。
#if defined(__linux__)
#  define CP_ALGO_USE_MMAP 1
#  include <sys/mman.h>
#else
#  define CP_ALGO_USE_MMAP 0
#endif

namespace cp_algo {
    template <typename T, size_t Align = 32>
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

    template<typename T> using big_vector = std::vector<T, big_alloc<T>>;
    template<typename T> using big_basic_string = std::basic_string<T, std::char_traits<T>, big_alloc<T>>;
    template<typename T> using big_deque = std::deque<T, big_alloc<T>>;
    template<typename T> using big_stack = std::stack<T, big_deque<T>>;
    template<typename T> using big_queue = std::queue<T, big_deque<T>>;
    template<typename T> using big_priority_queue = std::priority_queue<T, big_vector<T>>;
    template<typename T> using big_forward_list = std::forward_list<T, big_alloc<T>>;
    using big_string = big_basic_string<char>;

    template<typename Key, typename Value, typename Compare = std::less<Key>>
    using big_map = std::map<Key, Value, Compare, big_alloc<std::pair<const Key, Value>>>;
    template<typename T, typename Compare = std::less<T>>
    using big_multiset = std::multiset<T, Compare, big_alloc<T>>;
    template<typename T, typename Compare = std::less<T>>
    using big_set = std::set<T, Compare, big_alloc<T>>;
#ifdef CP_ALGO_HAS_GENERATOR
    template<typename Ref, typename V = void>
    using big_generator = std::generator<Ref, V, big_alloc<std::byte>>;
#endif
}

#ifdef CP_ALGO_HAS_GENERATOR
namespace std::ranges {
    template<typename Ref, typename V>
    elements_of(cp_algo::big_generator<Ref, V>&&) -> elements_of<cp_algo::big_generator<Ref, V>&&, cp_algo::big_alloc<std::byte>>;
}
#endif






#include <bit>

#if defined(__x86_64__) && !defined(CP_ALGO_DISABLE_AVX2) && !defined(USE_SIMDE)
#define CP_ALGO_BIT_OPS_TARGET _Pragma("GCC target(\"avx2,bmi,bmi2,lzcnt,popcnt\")")
#elif defined(__x86_64__)
#define CP_ALGO_BIT_OPS_TARGET _Pragma("GCC target(\"bmi,bmi2,lzcnt,popcnt\")")
#else
#define CP_ALGO_BIT_OPS_TARGET   /* non-x86: no target pragma */
#endif

#define CP_ALGO_BIT_PRAGMA_PUSH \
    _Pragma("GCC push_options") \
    CP_ALGO_BIT_OPS_TARGET

CP_ALGO_BIT_PRAGMA_PUSH
namespace cp_algo {
    template<typename Uint>
    constexpr size_t bit_width = sizeof(Uint) * 8;

    // n < 64
    uint64_t mask(size_t n) {
        return (1ULL << n) - 1;
    }
    size_t order_of_bit(auto x, size_t k) {
        return k ? std::popcount(x << (bit_width<decltype(x)> - k)) : 0;
    }
#ifdef __BMI2__
    inline size_t kth_set_bit(uint64_t x, size_t k) {
        return std::countr_zero(_pdep_u64(1ULL << k, x));
    }
#endif
    template<int fl = 0>
    void with_bit_floor(size_t n, auto &&callback) {
        if constexpr (fl >= 63) {
            return;
        } else if (n >> (fl + 1)) {
            with_bit_floor<fl + 1>(n, callback);
        } else {
            callback.template operator()<1ULL << fl>();
        }
    }
    void with_bit_ceil(size_t n, auto &&callback) {
        with_bit_floor(n, [&]<size_t N>() {
            if(N == n) {
                callback.template operator()<N>();
            } else {
                callback.template operator()<N << 1>();
            }
        });
    }

    // bitset 用の read_bits / write_bits 系は subset_exp 経路で未使用なので削除。
}
#pragma GCC pop_options





#include <chrono>

namespace cp_algo {
#ifdef CP_ALGO_CHECKPOINT
    big_map<big_string, double> checkpoints;
    double last;
#endif
    template<bool final = false>
    void checkpoint([[maybe_unused]] auto const& _msg) {
#ifdef CP_ALGO_CHECKPOINT
        big_string msg = _msg;
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
    template<bool final = false>
    void checkpoint() {
        checkpoint<final>("");
    }
}


#include <ranges>
#include <algorithm>

CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math {
    const size_t max_logn = 20;
    
    enum transform_dir { forw, inv };
    
    template<auto N, transform_dir direction>
    inline void xor_transform(auto &&a) {
        [[gnu::assume(N <= 1 << 30)]];
        if constexpr (N <= 32) {
            for (size_t i = 1; i < N; i *= 2) {
                for (size_t j = 0; j < N; j += 2 * i) {
                    for (size_t k = j; k < j + i; k++) {
                        for (size_t z = 0; z < max_logn; z++) {
                            auto x = a[k][z] + a[k + i][z];
                            auto y = a[k][z] - a[k + i][z];
                            a[k][z] = x;
                            a[k + i][z] = y;
                        }
                    }
                }
            }
        } else {
            auto add = [&](auto &a, auto &b) __attribute__((always_inline)) {
                auto x = a + b, y = a - b;
                a = x, b = y;
            };
            constexpr auto quar = N / 4;

            for (size_t i = 0; i < (size_t)quar; i++) {
                auto x0 = a[i + (size_t)quar * 0];
                auto x1 = a[i + (size_t)quar * 1];
                auto x2 = a[i + (size_t)quar * 2];
                auto x3 = a[i + (size_t)quar * 3];

                #pragma GCC unroll max_logn
                for (size_t z = 0; z < max_logn; z++) {
                    add(x0[z], x2[z]);
                    add(x1[z], x3[z]);
                }
                #pragma GCC unroll max_logn
                for (size_t z = 0; z < max_logn; z++) {
                    add(x0[z], x1[z]);
                    add(x2[z], x3[z]);
                }

                a[i + (size_t)quar * 0] = x0;
                a[i + (size_t)quar * 1] = x1;
                a[i + (size_t)quar * 2] = x2;
                a[i + (size_t)quar * 3] = x3;
            }
            xor_transform<quar, direction>(&a[quar * 0]);
            xor_transform<quar, direction>(&a[quar * 1]);
            xor_transform<quar, direction>(&a[quar * 2]);
            xor_transform<quar, direction>(&a[quar * 3]);
        }
    }

    
    template<transform_dir direction>
    inline void xor_transform(auto &&a, auto n) {
        with_bit_floor(n, [&]<auto NN>() {
            assert(NN == n);
            xor_transform<NN, direction>(a);
        });
    }
    
    template<transform_dir direction = forw>
    inline void xor_transform(auto &&a) {
        xor_transform<direction>(a, std::size(a));
    }

    // Generic rank vectors processor with variadic inputs
    // Assumes output[0] = 0, caller is responsible for handling rank 0
    // Returns the output array
    auto on_rank_vectors(auto &&cb, auto const& ...inputs) {
        static_assert(sizeof...(inputs) >= 1, "on_rank_vectors requires at least one input");
        
        // Create tuple of input references once
        auto input_tuple = std::forward_as_tuple(inputs...);
        auto const& first_input = std::get<0>(input_tuple);
        using base = std::decay_t<decltype(first_input[0])>;
        big_vector<base> out(std::size(first_input));
        
        auto N = std::size(first_input);
        constexpr size_t K = 4;
        N = std::max(N, 2 * K);
        const size_t n = std::bit_width(N) - 1;
        const size_t T = std::min<size_t>(n - 3, 2);
        const size_t bottoms = 1 << (n - T - 1);
        const auto M = std::size(first_input);
        
        // Create array buffers for each input
        auto create_buffers = [bottoms]<typename... Args>(const Args&...) {
            return std::make_tuple(
                big_vector<std::array<typename std::decay_t<Args>::value_type, max_logn>>(bottoms)...
            );
        };
        auto buffers = std::apply(create_buffers, input_tuple);
        
        checkpoint("alloc buffers");
        big_vector<uint32_t> counts(2 * bottoms);
        for(size_t i = 1; i < 2 * bottoms; i++) {
            counts[i] = (uint32_t)std::popcount(i);
        }
        checkpoint("prepare");
        
        for(size_t top = 0; top < N / 2; top += bottoms) {
            // Clear all buffers
            std::apply([bottoms](auto&... bufs) {
                (..., memset(bufs.data(), 0, sizeof(bufs[0]) * bottoms));
            }, buffers);
            checkpoint("memset");
            
            // Initialize buffers from inputs
            std::apply([&](auto const&... inps) {
                std::apply([&](auto&... bufs) {
                    auto init_one = [&](auto const& inp, auto& buf) {
                        for(size_t i = 0; i < M; i += 2 * bottoms) {
                            bool parity = __builtin_parity(uint32_t((i >> 1) & top));
                            size_t limit = std::min(M, i + 2 * bottoms) - i;
                            uint32_t count = (uint32_t)std::popcount(i) - 1;
                            for(size_t bottom = (i == 0); bottom < limit; bottom++) {
                                if (parity) {
                                    buf[bottom >> 1][count + counts[bottom]] -= inp[i + bottom];
                                } else {
                                    buf[bottom >> 1][count + counts[bottom]] += inp[i + bottom];
                                }
                            }
                        }
                    };
                    (init_one(inps, bufs), ...);
                }, buffers);
            }, input_tuple);
            
            checkpoint("init");
            std::apply([](auto&... bufs) {
                (..., xor_transform(bufs));
            }, buffers);
            checkpoint("transform");
            
            assert(bottoms % K == 0);
            for(size_t i = 0; i < bottoms; i += K) {
                std::apply([&](auto&... bufs) {
                    auto extract_one = [&](auto& buf) {
                        std::array<u64x4, max_logn> aa;
                        for(size_t j = 0; j < max_logn; j++) {
                            for(size_t z = 0; z < K; z++) {
                                aa[j][z] = buf[i + z][j].getr();
                            }
                        }
                        return aa;
                    };
                    
                    auto aa_tuple = std::make_tuple(extract_one(bufs)...);
                    std::apply(cb, aa_tuple);
                    
                    // Write results back: only first array needs to be written
                    auto& first_buf = std::get<0>(std::forward_as_tuple(bufs...));
                    const auto& first_aa = std::get<0>(aa_tuple);
                    for(size_t j = 0; j < max_logn; j++) {
                        for(size_t z = 0; z < K; z++) {
                            first_buf[i + z][j].setr((uint32_t)first_aa[j][z]);
                        }
                    }
                }, buffers);
            }
            
            checkpoint("dot");
            auto& first_buf = std::get<0>(buffers);
            xor_transform<inv>(first_buf);
            checkpoint("transform");
            
            // Gather results from first buffer

            for(size_t i = 0; i < M; i += 2 * bottoms) {
                bool parity = __builtin_parity(uint32_t((i >> 1) & top));
                size_t limit = std::min(M, i + 2 * bottoms) - i;
                uint32_t count = (uint32_t)std::popcount(i) - 1;
                for(size_t bottom = (i == 0); bottom < limit; bottom++) {
                    if (parity) {
                        out[i + bottom] -= first_buf[bottom >> 1][count + counts[bottom]];
                    } else {
                        out[i + bottom] += first_buf[bottom >> 1][count + counts[bottom]];
                    }
                }
            }
            checkpoint("gather");
        }
        const base ni = base(N / 2).inv();
        for(auto& x : out) {x *= ni;}
        return out;
    }

    template<typename base>
    big_vector<base> subset_convolution(std::span<base> f, std::span<base> g) {
        big_vector<base> outpa;
        with_bit_floor(std::size(f), [&]<auto N>() {
            constexpr size_t lgn = std::bit_width(N) - 1;
            [[gnu::assume(lgn <= max_logn)]];
            outpa = on_rank_vectors([](auto &a, auto const& b) {
                std::decay_t<decltype(a)> res = {};
                const auto mod = base::mod();
                const auto imod = math::inv2(-mod);
                const auto r4 = u64x4() + uint64_t(-1) % mod + 1;
                auto add = [&](size_t i) {
                    if constexpr (lgn) for(size_t j = 0; i + j + 1 < lgn; j++) {
                        res[i + j + 1] += (u64x4)_mm256_mul_epu32(__m256i(a[i]), __m256i(b[j]));
                    }
                };
                if constexpr (lgn) for(size_t i = 0; i < lgn; i++) { add(i); }
                if constexpr (lgn) if constexpr (lgn) for(size_t k = 0; k < lgn; k++) {
                    res[k] = montgomery_reduce(res[k], mod, imod);
                    res[k] = montgomery_mul(res[k], r4, mod, imod);
                    a[k] = res[k] >= mod ? res[k] - mod : res[k];
                }
            }, f, g);
            
            outpa[0] = f[0] * g[0];
            for(size_t i = 1; i < std::size(f); i++) {
                outpa[i] += f[i] * g[0] + f[0] * g[i];
            }
            checkpoint("fix 0");
        });
        return outpa;
    }

    template<typename base>
    big_vector<base> subset_exp(std::span<base> g) {
        if (size(g) == 1) {
            return big_vector<base>{1};
        }
        size_t N = std::size(g);
        auto out0 = subset_exp(std::span(g).first(N / 2));
        auto out1 = subset_convolution<base>(out0, std::span(g).last(N / 2));
        out0.insert(end(out0), begin(out1), end(out1));
        cp_algo::checkpoint("extend out");
        return out0;
    }

    // subset_compose / subset_conv_transpose / subset_power_projection は本問題で未使用、
    // libc++ の views::enumerate / 大域 CTAD で詰まるので削除。
}
#pragma GCC pop_options


// ---------------------------------------------------------------------------
// Wrapper
// ---------------------------------------------------------------------------
struct SubsetExp {
 static vector<u32> run(int N, const vector<u32>& b_in) {
  using base = cp_algo::math::modint<(int64_t) MOD>;
  size_t sz = (size_t) 1 << N;
  cp_algo::big_vector<base> a(sz);
  for (size_t i = 0; i < sz; ++i) a[i].setr((typename base::UInt) b_in[i]);
  auto c = cp_algo::math::subset_exp<base>(std::span<base>(a));
  vector<u32> result(sz);
  for (size_t i = 0; i < sz; ++i) {
   u32 v = (u32) c[i].getr();
   if (v >= MOD) v -= MOD;
   result[i] = v;
  }
  return result;
 }
};
