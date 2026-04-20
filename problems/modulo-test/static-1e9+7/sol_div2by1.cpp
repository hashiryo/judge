// https://gmplib.org/~tege/division-paper.pdf
#include <bits/stdc++.h>
using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;
struct MP_D2B1 {  // mod < 2^31
 u32 mod;
 constexpr MP_D2B1(): mod(0), s(0), d(0), v(0) {}
 constexpr MP_D2B1(u32 m): mod(m), s(__builtin_clz(m)), d(m << s), v(u64(-1) / d) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 constexpr inline u32 set(u32 n) const { return n; }
 constexpr inline u32 get(u32 n) const { return n; }
 constexpr inline u32 norm(u32 n) const { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const { return l+= r, l < mod ? l : l - mod; }
 constexpr inline u32 diff(u32 l, u32 r) const { return l-= r, l >> 31 ? l + mod : l; }
private:
 u8 s;
 u32 d;
 u32 v;
 constexpr inline u32 rem(u64 n) const {
  n<<= s;
  u64 q= u64(v) * u32(n >> 32) + n;
  u32 q1= u32(q >> 32) + 1;
  u32 r= n - q1 * d;
  if(r > u32(q)) r+= d;
  if(r >= d) r-= d;
  return r >> s;
 }
};
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 constexpr MP_D2B1 mp(int(1e9 + 7));
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
