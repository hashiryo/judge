#include <bits/stdc++.h>
using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;
struct MP_Br {  // mod < 2^31
 u32 mod;
 constexpr MP_Br(): mod(0), s(0), x(0) {}
 constexpr MP_Br(u32 m): mod(m), s(95 - __builtin_clz(m - 1)), x(((u128(1) << s) - 1) / m + 1) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 static constexpr inline u32 set(u32 n) { return n; }
 static constexpr inline u32 get(u32 n) { return n; }
 static constexpr inline u32 norm(u32 n) { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> 31 ? l + mod : l; }
private:
 u8 s;
 u64 x;
 constexpr inline u32 quo(u64 n) const { return (u128(x) * n) >> s; }
 constexpr inline u32 rem(u64 n) const { return n - u64(quo(n)) * mod; }
};
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 constexpr MP_Br mp((1u << 31) - 1);
 u32 n, state, a, b, c, d;
 cin >> n >> state >> a >> b >> c >> d;
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
