// https://zenn.dev/yatyou/articles/f8a3969c6e9f7b
#include <bits/stdc++.h>
using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;
struct MP_Pl {  // mod < 2^32/phi
 u32 mod;
 constexpr MP_Pl(): mod(0), iv(0), r2(0) {}
 constexpr MP_Pl(u32 m): mod(m), iv(inv(m)), r2(-u128(mod) % mod) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return reduce(u64(l) * r); }
 constexpr inline u32 set(u32 n) const { return mul(n, r2); }
 constexpr inline u32 get(u32 n) const { return reduce(n); }
 constexpr inline u32 norm(u32 n) const { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> 31 ? l + mod : l; }
private:
 u64 iv, r2;
 static constexpr u64 inv(u64 n, int e= 6, u64 x= 1) { return e ? inv(n, e - 1, x * (2 - x * n)) : x; }
 constexpr inline u32 reduce(u64 w) const { return (u128((w * iv) | u32(-1)) * mod) >> 64; }
};
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 u32 n, state, a, b, c, d;
 cin >> n >> state >> a >> b >> c >> d;
 MP_Pl mp(998244353);
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
