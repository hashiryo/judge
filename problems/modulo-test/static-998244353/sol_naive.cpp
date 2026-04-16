#include <bits/stdc++.h>
using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;
struct MP_Na {  // mod < 2^32
 u32 mod;
 constexpr MP_Na(): mod(0) {}
 constexpr MP_Na(u32 m): mod(m) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return u64(l) * r % mod; }
 constexpr inline u32 set(u32 n) const { return n; }
 constexpr inline u32 get(u32 n) const { return n; }
 constexpr inline u32 norm(u32 n) const { return n; }
 constexpr inline u32 plus(u64 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u64 l, u32 r) const { return l-= r, l >> 63 ? l + mod : l; }
};
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 constexpr MP_Na mp(998244353);
 u32 n, state, a, b;
 cin >> n >> state >> a >> b;
 state= mp.set(state);
 a= mp.set(a);
 b= mp.set(b);
 for(int i= 0; i < n; ++i) state= mp.plus(mp.mul(state, a), b);
 state= mp.get(state);
 cout << state << '\n';
 return 0;
}
