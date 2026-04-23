#include <bits/stdc++.h>
using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;
struct MP_Br {  // mod < 2^30
 u32 mod;
 MP_Br(): mod(0), x(0) {}
 MP_Br(u32 m): mod(m), x(-u64(m) / m + 1) { asm volatile("" : "+r"(mod), "+r"(x)); }
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 static constexpr inline u32 set(u32 n) { return n; }
 constexpr inline u32 get(u32 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u32 norm(u32 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < (mod << 1) ? l : l - (mod << 1); }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> 31 ? l + (mod << 1) : l; }
private:
 u64 x;
 constexpr inline u32 quo(u64 n) const { return (u128(n) * x) >> 64; }
 constexpr inline u32 rem(u64 n) const { return n - u64(quo(n)) * mod; }
};
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 u32 n, state, a, b, c, d;
 cin >> n >> state >> a >> b >> c >> d;
 MP_Br mp(998244353);
 state= mp.set(state);
 a= mp.set(a);
 b= mp.set(b);
 c= mp.set(c);
 d= mp.set(d);
 for(int i= 0; i < n; ++i) state= mp.mul(mp.plus(mp.mul(state, a), b), mp.plus(mp.mul(state, c), d));
 state= mp.get(state);
 cout << state << '\n';
 return 0;
}
