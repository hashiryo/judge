// https://gmplib.org/~tege/division-paper.pdf
#include <bits/stdc++.h>
using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;
struct MP_D2B1 {  // 2^31 <= mod < 2^32
 u32 mod;
 constexpr MP_D2B1(): mod(0), v(0) {}
 constexpr MP_D2B1(u32 m): mod(m), v(u64(-1) / m) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 constexpr inline u32 set(u32 n) const { return n; }
 constexpr inline u32 get(u32 n) const { return n; }
 constexpr inline u32 norm(u32 n) const { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const {
  l+= r;
  if(l < r) l-= mod;
  if(l >= mod) l-= mod;
  return l;
 }
 constexpr inline u32 diff(u32 l, u32 r) const { return plus(l, mod - r); }
private:
 u32 v;
 constexpr inline u32 rem(u64 n) const {
  u64 q= u64(v) * u32(n >> 32) + n;
  u32 q1= u32(q >> 32) + 1;
  u32 r= n - q1 * mod;
  r+= mod & -u32(r > u32(q));
  if(__builtin_expect(r >= mod, 0)) r-= mod;
  return r;
 }
};
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 u32 n, mod, state, a, b, c, d;
 cin >> n >> mod >> state >> a >> b >> c >> d;
 MP_D2B1 mp(mod);
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
