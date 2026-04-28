#pragma once
#include "_common.hpp"
// https://gmplib.org/~tege/divcnst-pldi94.pdf section 7
using f80= long double;
struct MP {  // mod < 2^{30.5}
 u32 mod;
 constexpr MP(): mod(0) {}
 constexpr MP(u32 m): mod(m) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 static constexpr inline u32 set(u32 n) { return n; }
 static constexpr inline u32 get(u32 n) { return n; }
 static constexpr inline u32 norm(u32 n) { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> 31 ? l + mod : l; }
private:
 static constexpr f80 EPS= 1.0L + 0x1.0p-62L;
 constexpr inline u64 quo(u64 n) const { return EPS * n / mod; }
 constexpr inline u32 rem(u64 n) const { return n - quo(n) * mod; }
};
