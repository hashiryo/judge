// https://gmplib.org/~tege/divcnst-pldi94.pdf section 7
#include <bits/stdc++.h>
using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;
using f80= long double;
struct MP_F80 {  // mod < 2^{30.5}
 u32 mod;
 constexpr MP_F80(): mod(0) {}
 constexpr MP_F80(u32 m): mod(m) {}
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
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 u32 n, mod, state, a, b;
 cin >> n >> mod >> state >> a >> b;
 MP_F80 mp(mod);
 state= mp.set(state);
 a= mp.set(a);
 b= mp.set(b);
 for(int i= 0; i < n; ++i) state= mp.plus(mp.mul(state, a), b);
 state= mp.get(state);
 cout << state << '\n';
 return 0;
}
