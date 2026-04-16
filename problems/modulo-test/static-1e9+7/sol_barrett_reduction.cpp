#include <bits/stdc++.h>
using namespace std;
using u8= unsigned char;
using u32= unsigned;
using i64= long long;
using u64= unsigned long long;
using u128= __uint128_t;
struct MP_Br {  // 2^20 < mod <= 2^41
 u64 mod;
 constexpr MP_Br(): mod(0), x(0) {}
 constexpr MP_Br(u64 m): mod(m), x((u128(1) << 84) / m) {}
 constexpr inline u64 mul(u64 l, u64 r) const { return rem(u128(l) * r); }
 static constexpr inline u64 set(u64 n) { return n; }
 constexpr inline u64 get(u64 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u64 norm(u64 n) const { return n >= mod ? n - mod : n; }
 constexpr inline u64 plus(u64 l, u64 r) const { return l+= r, l < (mod << 1) ? l : l - (mod << 1); }
 constexpr inline u64 diff(u64 l, u64 r) const { return l-= r, l >> 63 ? l + (mod << 1) : l; }
private:
 u64 x;
 constexpr inline u128 quo(const u128& n) const { return (n * x) >> 84; }
 constexpr inline u64 rem(const u128& n) const { return n - quo(n) * mod; }
};
signed main() {
 cin.tie(0);
 ios::sync_with_stdio(false);
 constexpr MP_Br mp(int(1e9 + 7));
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
