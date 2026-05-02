#pragma once
// =============================================================================
// algos/bpsw.hpp の判定ロジックを保ったまま、入力 4 個ぶんの Lucas chain を
// インターリーブして進める "batch=4" 版。
//
// 動機:
//   bpsw.hpp の Lucas 主ループは um/vm/qn の RAW 依存で直線律速 (IPC ~1.0-1.5)。
//   1 イテレーション内の独立 mrmul を 4 並列にし、CPU の mul port を埋める
//   ことで wallclock を縮める。判定ロジック自体は同じ BPSW (MR-2 + 強 Lucas)。
//
// 戦略:
//   - Mont/Jacobi/MR-2/Lucas 主体は bpsw.hpp と同一実装
//   - 入力を B=4 個ずつ chunk
//   - 各 chunk で per-slot 初期化 (D 探索, dm/qm 計算) を実行
//   - Lucas 主ループは "while any k[i]" で全 slot 同時に doubling+1step
//     active な slot だけ state を update (compiler が cmov 化)
//   - Strong suffix も同様に max iter で全 slot 同時
//
// 注意:
//   - lucas_dprobe (D 探索) は per-slot で逐次 (Jacobi 計算は分岐数が
//     入力依存なので並列化の旨味が薄い + D が決まらないと chain が回せない)
//   - small case (n<4 / 偶数 / 小素数で除算) は batch に入れず即決
// =============================================================================

#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"

namespace yosupo_bpsw_batch4_mr2 {
using std::uint32_t;
using std::uint64_t;
using std::int64_t;
using u128_local = __uint128_t;

// =============================================================================
// Mont, Jacobi: bpsw.hpp と同一 (コピー)
// =============================================================================
struct Mont {
 uint64_t m;
 uint64_t mi;
 uint64_t mh;
 uint64_t R;
 uint64_t Rn;
 uint64_t R2;
 uint64_t d;
 int k;
 Mont() = default;
 explicit Mont(uint64_t mm): m(mm) {
  uint64_t x = mm * 3 ^ 2;
  uint64_t y = (uint64_t) 1 - mm * x;
  for (int w = 5; w < 64; w *= 2) { x = x * (y + 1); y = y * y; }
  mi = x;
  mh = (mm >> 1) + 1;
  R = (uint64_t) -(int64_t) 1 % mm + 1;
  R2 = (uint64_t) ((u128_local) -(__int128) 1 % mm) + 1;
  Rn = mm - R;
  d = mm - 1;
  k = __builtin_ctzll(d);
  d >>= k;
 }
 inline uint64_t add(uint64_t a, uint64_t b) const {
  uint64_t s = a + b;
  return s >= m ? s - m : s;
 }
 inline uint64_t sub(uint64_t a, uint64_t b) const {
  return a >= b ? a - b : a + (m - b);
 }
 inline uint64_t div2(uint64_t ar) const {
  return (ar >> 1) + ((ar & 1) ? mh : 0);
 }
 inline uint64_t mrmul(uint64_t ar, uint64_t br) const {
  u128_local t = (u128_local) ar * br;
  uint64_t q = (uint64_t) t * mi;
  uint64_t hi = (uint64_t) (t >> 64);
  uint64_t mq_hi = (uint64_t) (((u128_local) q * m) >> 64);
  uint64_t res = hi - mq_hi;
  return hi < mq_hi ? res + m : res;
 }
 inline uint64_t to_mont(uint64_t a) const { return mrmul(a, R2); }
 inline uint64_t mod_n(uint64_t x) const { return x < m ? x : x % m; }
 inline uint64_t pow_mont(uint64_t ar, uint64_t b) const {
  if (b == 0) return R;
  while ((b & 1) == 0) { ar = mrmul(ar, ar); b >>= 1; }
  uint64_t t = ar; b >>= 1;
  while (b) {
   ar = mrmul(ar, ar);
   if (b & 1) t = mrmul(t, ar);
   b >>= 1;
  }
  return t;
 }
};

constexpr int B = 4;

// MR-2 scalar (Pass 1 → Pass 2 で 1 入力だけのときの fallback 用)
inline bool mr_base2(const Mont& mo) {
 uint64_t br = mo.add(mo.R, mo.R);
 br = mo.pow_mont(br, mo.d);
 if (br == mo.R || br == mo.Rn) return true;
 for (int i = 1; i < mo.k; ++i) {
  br = mo.mrmul(br, br);
  if (br == mo.Rn) return true;
 }
 return false;
}

// =============================================================================
// MR-2 batch=4。
// 4 入力分の base-2 強擬素数判定を並列に走らせる。
// pow_mont は (br = 2^d mod m) を求めるが、これを 4 slot の squaring を
// インターリーブして実行 → mul ユニットを埋める。
//
// 出力: out[i] が true なら "base-2 で probable prime" (Lucas に進める)。
//       false なら確定 composite。
// active[i] が false の slot は触らない。
// =============================================================================
__attribute__((always_inline))
inline void mr2_pow_step_one(
    const Mont& mo, uint64_t& ar_io, uint64_t& t_io,
    uint64_t& b_io, bool& seeded_io)
{
 uint64_t b = b_io;
 if (b == 0) return;
 if (!seeded_io) {
  // 「最下位 1 が見つかるまで ar を squaring」フェーズ
  if ((b & 1) == 0) {
   ar_io = mo.mrmul(ar_io, ar_io);
   b_io = b >> 1;
   return;
  }
  // 最初の 1 bit に当たった: t = ar, b >>= 1
  t_io = ar_io;
  b_io = b >> 1;
  seeded_io = true;
  return;
 }
 // 通常フェーズ: 毎反復 ar squaring、bit が立ってたら t *= ar
 ar_io = mo.mrmul(ar_io, ar_io);
 if (b & 1) t_io = mo.mrmul(t_io, ar_io);
 b_io = b >> 1;
}

inline void mr_base2_batch(const Mont mos[B], bool active[B], bool out[B]) {
 // pow_mont の状態: ar (Montgomery 形式の base 2)、t (累積結果)、b (残り指数)、seeded
 uint64_t ar0=0,ar1=0,ar2=0,ar3=0;
 uint64_t t0=0,t1=0,t2=0,t3=0;
 uint64_t b0=0,b1=0,b2=0,b3=0;
 bool seeded0=false, seeded1=false, seeded2=false, seeded3=false;

 auto init = [](const Mont& mo, bool act, uint64_t& ar, uint64_t& t,
                uint64_t& b, bool& seeded) {
  if (!act) return;
  ar = mo.add(mo.R, mo.R); // base 2 in Montgomery
  t = mo.R;                 // 初期値 (使用前に t = ar で上書き)
  b = mo.d;
  seeded = false;
  if (b == 0) {
   // 通常はあり得ないが念のため (b == 0 なら結果 = R)
   t = mo.R;
   seeded = true;
  }
 };
 init(mos[0], active[0], ar0, t0, b0, seeded0);
 init(mos[1], active[1], ar1, t1, b1, seeded1);
 init(mos[2], active[2], ar2, t2, b2, seeded2);
 init(mos[3], active[3], ar3, t3, b3, seeded3);

 // 4 slot の pow_mont をインターリーブ
 while ((b0 | b1 | b2 | b3) != 0) {
  mr2_pow_step_one(mos[0], ar0, t0, b0, seeded0);
  mr2_pow_step_one(mos[1], ar1, t1, b1, seeded1);
  mr2_pow_step_one(mos[2], ar2, t2, b2, seeded2);
  mr2_pow_step_one(mos[3], ar3, t3, b3, seeded3);
 }
 // この時点で t[i] = 2^d mod m (Montgomery 形式)
 uint64_t br0 = t0, br1 = t1, br2 = t2, br3 = t3;

 // 強判定: br == R or br == Rn なら prime probable。
 // それ以外で br を k-1 回 squaring して、その途中で Rn に至れば prime。
 // 各 slot の k が違うので per-slot に reservoir 持つ。
 int rem0 = active[0] ? mos[0].k - 1 : 0;
 int rem1 = active[1] ? mos[1].k - 1 : 0;
 int rem2 = active[2] ? mos[2].k - 1 : 0;
 int rem3 = active[3] ? mos[3].k - 1 : 0;
 bool decided0=!active[0], decided1=!active[1], decided2=!active[2], decided3=!active[3];

 auto first_check = [](const Mont& mo, bool act, uint64_t br,
                       bool& decided, bool& out_i, int& rem) {
  if (!act) return;
  if (br == mo.R || br == mo.Rn) {
   out_i = true;
   decided = true;
   rem = 0; // 強判定 squaring は不要
  }
 };
 first_check(mos[0], active[0], br0, decided0, out[0], rem0);
 first_check(mos[1], active[1], br1, decided1, out[1], rem1);
 first_check(mos[2], active[2], br2, decided2, out[2], rem2);
 first_check(mos[3], active[3], br3, decided3, out[3], rem3);

 // 残ったやつだけ squaring を per-slot reservoir で並列に進める。
 // 通常 k は 1〜数 bit (ctz(m-1)) なので回数は少ないが、batch でやると
 // pass 2 全体の独立 mul を稼ぐのに地味に効く。
 while (rem0 + rem1 + rem2 + rem3 > 0) {
  if (!decided0 && rem0 > 0) {
   br0 = mos[0].mrmul(br0, br0);
   if (br0 == mos[0].Rn) { out[0] = true; decided0 = true; rem0 = 0; }
   else rem0--;
  }
  if (!decided1 && rem1 > 0) {
   br1 = mos[1].mrmul(br1, br1);
   if (br1 == mos[1].Rn) { out[1] = true; decided1 = true; rem1 = 0; }
   else rem1--;
  }
  if (!decided2 && rem2 > 0) {
   br2 = mos[2].mrmul(br2, br2);
   if (br2 == mos[2].Rn) { out[2] = true; decided2 = true; rem2 = 0; }
   else rem2--;
  }
  if (!decided3 && rem3 > 0) {
   br3 = mos[3].mrmul(br3, br3);
   if (br3 == mos[3].Rn) { out[3] = true; decided3 = true; rem3 = 0; }
   else rem3--;
  }
 }
 if (active[0] && !decided0) out[0] = false;
 if (active[1] && !decided1) out[1] = false;
 if (active[2] && !decided2) out[2] = false;
 if (active[3] && !decided3) out[3] = false;
}

// Jacobi (bpsw.hpp と同一)
inline int jacobi_u64(uint64_t a, uint64_t n) {
 int j = 1;
 while (a) {
  int b = __builtin_ctzll(a);
  if ((b & 1) && ((n + 2) & 5) == 5) j = -j;
  a >>= b;
  if (a < n) {
   if (a > (n >> 1)) {
    a = n - a;
    if ((n & 3) == 3) j = -j;
    int b2 = __builtin_ctzll(a);
    if ((b2 & 1) && ((n + 2) & 5) == 5) j = -j;
    a >>= b2;
   }
   if ((a & n & 3) == 3) j = -j;
   uint64_t t = n - a;
   n = a;
   a = t;
  } else {
   a -= n;
  }
 }
 return n == 1 ? j : 0;
}
inline int jacobi_signed(int64_t ai, uint64_t n) {
 if (ai >= 0) return jacobi_u64((uint64_t) ai, n);
 int sign = ((n & 3) == 3) ? -1 : 1;
 return sign * jacobi_u64((uint64_t) -ai, n);
}

inline bool is_perfect_square(uint64_t n) {
 if (n == 0) return true;
 uint64_t s = (uint64_t) std::sqrt((double) n);
 for (uint64_t t = (s ? s - 1 : 0); t <= s + 1; ++t) if (t * t == n) return true;
 return false;
}

// D = 5, -7, 9, -11, ...; jacobi(D, n) == -1 の最初の D を返す。0 → 合成数。
inline int lucas_dprobe(uint64_t n) {
 for (uint32_t mag = 5;; mag += 2) {
  int sign = (mag & 2) ? -1 : 1;
  int64_t D = (int64_t) mag * sign;
  int j = jacobi_signed(D, n);
  if (j == -1) return (int) D;
  if (j == 0 && mag < n) return 0;
  if (mag == 61) {
   if (is_perfect_square(n)) return 0;
  }
  if (mag > 200000) return 0;
 }
}

// =============================================================================
// Batch=4 の Lucas chain。
//
// 入力:
//   mos[i]:    各 slot の Montgomery context
//   active[i]: false なら slot i は処理せず out[i] を保つ (= 既に決着済み)
// 出力:
//   out[i]:    active=true の slot のみ更新 (true=probable prime, false=composite)
//
// 主ループは 4 slot 分の doubling+1step を毎反復計算し、active && (k!=0)
// な slot だけ state を確定する。compiler は active マスクを cmov 化して
// 4 並列に独立な mrmul を発行するため ILP が稼げる。
// =============================================================================

// 1 slot ぶんの doubling + (top bit set 時) 1-step を行う inline ヘルパ。
// 各 slot は独立な local 変数を扱うので、4 回 inline 展開すれば
// コンパイラは 4 並列の mrmul を独立な命令列として並べられる。
__attribute__((always_inline))
inline void lucas_step_one(
    const Mont& mo, uint64_t qm_i, uint64_t dm_i,
    uint64_t& um_io, uint64_t& vm_io, uint64_t& qn_io, uint64_t& k_io)
{
 uint64_t k_i = k_io;
 if (k_i == 0) return;
 uint64_t um_i = um_io, vm_i = vm_io, qn_i = qn_io;
 // doubling
 uint64_t um_d = mo.mrmul(um_i, vm_i);
 uint64_t vm_d = mo.sub(mo.mrmul(vm_i, vm_i), mo.add(qn_i, qn_i));
 uint64_t qn_d = mo.mrmul(qn_i, qn_i);
 // 1-step は条件付き (bit が立ってるときだけ追加 mul)。eager 計算は
 // mul throughput を浪費するため、Apple M2 等の mul 帯域に余裕の少ない
 // ターゲットでは branch のほうが速い。bit は ~50% でしか立たないので
 // 平均 mul 数を ~33% 削減できる。
 if ((int64_t) k_i < 0) {
  uint64_t uu = mo.div2(mo.add(um_d, vm_d));
  uint64_t vv = mo.div2(mo.add(mo.mrmul(dm_i, um_d), vm_d));
  um_io = uu;
  vm_io = vv;
  qn_io = mo.mrmul(qn_d, qm_i);
 } else {
  um_io = um_d;
  vm_io = vm_d;
  qn_io = qn_d;
 }
 k_io = k_i << 1;
}

__attribute__((always_inline))
inline void lucas_strong_one(
    const Mont& mo, uint64_t& um_io, uint64_t& vm_io, uint64_t& qn_io,
    uint64_t& x_io, bool& done_io, bool& out_io)
{
 uint64_t x = x_io;
 if (x == 0 || done_io) return;
 uint64_t um_i = um_io, vm_i = vm_io, qn_i = qn_io;
 uint64_t um_n = mo.mrmul(um_i, vm_i);
 uint64_t vm_n = mo.sub(mo.mrmul(vm_i, vm_i), mo.add(qn_i, qn_i));
 uint64_t qn_n = mo.mrmul(qn_i, qn_i);
 um_io = um_n;
 vm_io = vm_n;
 qn_io = qn_n;
 if (vm_n == 0) {
  out_io = true;
  done_io = true;
  x_io = 0;
  return;
 }
 x_io = x >> 1;
}

// Mont を value で持って alias 推論を簡単にする。null slot 用に dummy を使う。
inline void lucas_batch(const Mont mos_in[B], bool active[B], bool out[B]) {
 // 値コピーで持つ (alias を切る + cache 局所性)
 Mont mo0 = active[0] ? mos_in[0] : Mont(3); // dummy: 任意の奇数
 Mont mo1 = active[1] ? mos_in[1] : Mont(3);
 Mont mo2 = active[2] ? mos_in[2] : Mont(3);
 Mont mo3 = active[3] ? mos_in[3] : Mont(3);

 uint64_t um0=0,vm0=0,qn0=0,qm0=0,dm0=0,k0=0;
 uint64_t um1=0,vm1=0,qn1=0,qm1=0,dm1=0,k1=0;
 uint64_t um2=0,vm2=0,qn2=0,qm2=0,dm2=0,k2=0;
 uint64_t um3=0,vm3=0,qn3=0,qm3=0,dm3=0,k3=0;

 auto init_slot = [](const Mont& mo, bool& act, bool& out_i,
                     uint64_t& um, uint64_t& vm, uint64_t& qn,
                     uint64_t& qm, uint64_t& dm, uint64_t& k) {
  if (!act) return;
  uint64_t n = mo.m;
  int D = lucas_dprobe(n);
  if (D == 0) { out_i = false; act = false; return; }
  uint64_t Q = (D < 0) ? (uint64_t)(1 - D) / 4
                       : n - (uint64_t)(D - 1) / 4;
  qm = mo.to_mont(mo.mod_n(Q));
  dm = (D < 0) ? mo.to_mont(n - mo.mod_n((uint64_t) -D))
               : mo.to_mont(mo.mod_n((uint64_t) D));
  um = mo.R;
  vm = mo.R;
  qn = qm;
  k = (n + 1) << __builtin_clzll(n + 1);
  k <<= 1;
 };
 init_slot(mo0, active[0], out[0], um0, vm0, qn0, qm0, dm0, k0);
 init_slot(mo1, active[1], out[1], um1, vm1, qn1, qm1, dm1, k1);
 init_slot(mo2, active[2], out[2], um2, vm2, qn2, qm2, dm2, k2);
 init_slot(mo3, active[3], out[3], um3, vm3, qn3, qm3, dm3, k3);

 // 主ループ: 4 slot を完全に独立な命令列として展開。
 while ((k0 | k1 | k2 | k3) != 0) {
  lucas_step_one(mo0, qm0, dm0, um0, vm0, qn0, k0);
  lucas_step_one(mo1, qm1, dm1, um1, vm1, qn1, k1);
  lucas_step_one(mo2, qm2, dm2, um2, vm2, qn2, k2);
  lucas_step_one(mo3, qm3, dm3, um3, vm3, qn3, k3);
 }

 // 主ループ後: um==0 or vm==0 で確定 prime
 bool done0 = false, done1 = false, done2 = false, done3 = false;
 if (active[0] && (um0 == 0 || vm0 == 0)) { out[0] = true; done0 = true; }
 if (active[1] && (um1 == 0 || vm1 == 0)) { out[1] = true; done1 = true; }
 if (active[2] && (um2 == 0 || vm2 == 0)) { out[2] = true; done2 = true; }
 if (active[3] && (um3 == 0 || vm3 == 0)) { out[3] = true; done3 = true; }

 uint64_t x0 = 0, x1 = 0, x2 = 0, x3 = 0;
 auto init_strong = [](const Mont& mo, bool act, bool done, bool& out_i, uint64_t& x) {
  if (!act || done) return;
  uint64_t n = mo.m;
  uint64_t xv = (n + 1) & ~n;
  x = xv >> 1;
  out_i = false;
 };
 init_strong(mo0, active[0], done0, out[0], x0);
 init_strong(mo1, active[1], done1, out[1], x1);
 init_strong(mo2, active[2], done2, out[2], x2);
 init_strong(mo3, active[3], done3, out[3], x3);

 while ((x0 | x1 | x2 | x3) != 0) {
  lucas_strong_one(mo0, um0, vm0, qn0, x0, done0, out[0]);
  lucas_strong_one(mo1, um1, vm1, qn1, x1, done1, out[1]);
  lucas_strong_one(mo2, um2, vm2, qn2, x2, done2, out[2]);
  lucas_strong_one(mo3, um3, vm3, qn3, x3, done3, out[3]);
 }
 if (active[0] && !done0) out[0] = false;
 if (active[1] && !done1) out[1] = false;
 if (active[2] && !done2) out[2] = false;
 if (active[3] && !done3) out[3] = false;
}

// =============================================================================
// 1 入力ぶんの fast-path 判定: small / even / 小素数除算で決着するケース。
// 戻り値:
//   first  = 決着したか (true: out が確定値)
//   second = 決着時の判定 (true=prime, false=composite)
// =============================================================================
inline std::pair<bool, bool> fast_decide(uint64_t n) {
 if (n < 4) return {true, n >= 2};
 if ((n & 1) == 0) return {true, false};
 for (uint64_t p : {3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull, 23ull, 29ull, 31ull, 37ull, 41ull, 43ull}) {
  if (n % p == 0) return {true, n == p};
 }
 return {false, false};
}

} // namespace yosupo_bpsw_batch4_mr2

// 3-pass 構成:
//   pass 1: trial division (fast_decide) で 80%+ を決着
//   pass 2: MR-2 を per-slot で実行 (composite を早期に弾く)
//   pass 3: 生き残ったものだけを 4 個ずつ Lucas batch
//
// pass 3 で「常に 4 slot が active」な状態を作ることで、Lucas batch の
// 並列利得を fully amortize できる。雑にチャンクすると 86% の fast-path
// 入力で batch slot が無駄になり、scalar より遅くなる。
struct Primality {
 static vector<bool> run(const vector<u64>& qs) {
  using namespace yosupo_bpsw_batch4_mr2;
  const size_t Q = qs.size();
  vector<bool> ans(Q, false);

  // pass 1: trial division。生き残った入力は MR-2 batch 用に集約。
  std::vector<uint32_t> mr2_idx;
  std::vector<Mont> mr2_mos;
  mr2_idx.reserve(Q / 4);
  mr2_mos.reserve(Q / 4);
  for (size_t i = 0; i < Q; ++i) {
   auto [decided, v] = fast_decide(qs[i]);
   if (decided) {
    ans[i] = v;
   } else {
    mr2_idx.push_back((uint32_t) i);
    mr2_mos.push_back(Mont(qs[i]));
   }
  }

  // pass 2: MR-2 を 4 並列で。生き残ったものを Lucas batch 用に集約。
  std::vector<uint32_t> need_lucas_idx;
  std::vector<Mont> need_lucas_mos;
  need_lucas_idx.reserve(mr2_idx.size() / 4);
  need_lucas_mos.reserve(mr2_idx.size() / 4);
  {
   const size_t M = mr2_idx.size();
   size_t k = 0;
   for (; k + B <= M; k += B) {
    Mont mos[B] = {mr2_mos[k], mr2_mos[k+1], mr2_mos[k+2], mr2_mos[k+3]};
    bool active[B] = {true, true, true, true};
    bool out[B] = {false, false, false, false};
    mr_base2_batch(mos, active, out);
    for (int j = 0; j < B; ++j) {
     if (out[j]) {
      need_lucas_idx.push_back(mr2_idx[k + j]);
      need_lucas_mos.push_back(mos[j]);
     } else {
      ans[mr2_idx[k + j]] = false;
     }
    }
   }
   // tail: scalar mr_base2 で
   for (; k < M; ++k) {
    if (mr_base2(mr2_mos[k])) {
     need_lucas_idx.push_back(mr2_idx[k]);
     need_lucas_mos.push_back(mr2_mos[k]);
    } else {
     ans[mr2_idx[k]] = false;
    }
   }
  }

  // pass 3: Lucas batch 4 並列
  const size_t L = need_lucas_idx.size();
  size_t k = 0;
  for (; k + B <= L; k += B) {
   Mont mos[B] = {need_lucas_mos[k], need_lucas_mos[k+1],
                  need_lucas_mos[k+2], need_lucas_mos[k+3]};
   bool active[B] = {true, true, true, true};
   bool out[B] = {false, false, false, false};
   lucas_batch(mos, active, out);
   for (int j = 0; j < B; ++j) ans[need_lucas_idx[k + j]] = out[j];
  }
  // 末尾 (< B 個): scalar batch (slot 0 のみ active)
  for (; k < L; ++k) {
   Mont mos[B] = {need_lucas_mos[k], need_lucas_mos[k],
                  need_lucas_mos[k], need_lucas_mos[k]};
   bool active[B] = {true, false, false, false};
   bool out[B] = {false, false, false, false};
   lucas_batch(mos, active, out);
   ans[need_lucas_idx[k]] = out[0];
  }
  return ans;
 }
};
