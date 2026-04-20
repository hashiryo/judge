#include <bits/stdc++.h>
using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;
struct MP_Br {  // mod < 2^32
 u32 mod;
 constexpr MP_Br(): mod(0), x(0) {}
 constexpr MP_Br(u32 m): mod(m), x(((u128(1) << 96) - 1) / m + 1) {}
 constexpr inline u32 mul(u32 l, u32 r) const { return rem(u64(l) * r); }
 static constexpr inline u32 set(u32 n) { return n; }
 static constexpr inline u32 get(u32 n) { return n; }
 static constexpr inline u32 norm(u32 n) { return n; }
 constexpr inline u32 plus(u32 l, u32 r) const {
  l+= r;
  if(l < r) l-= mod;
  if(l >= mod) l-= mod;
  return l;
 }
 constexpr inline u32 diff(u32 l, u32 r) const { return plus(l, mod - r); }
private:
 u128 x;
 constexpr inline u32 quo(u64 n) const { return (x * n) >> 96; }
 constexpr inline u32 rem(u64 n) const { return n - u64(quo(n)) * mod; }
};
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 u32 n, mod, state, a, b, c, d;
 cin >> n >> mod >> state >> a >> b >> c >> d;
 MP_Br mp(mod);
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
