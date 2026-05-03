#pragma once
#include "_common.hpp"
// yosupo kth_root_integer 最速提出 (yzlf, https://judge.yosupo.jp/submission/189811)
// のアルゴリズム部のみを移植 (独自 I/O 部は省略、 ostringstream 経由)。
//
// アイデア:
//   k で分岐:
//     k = 1: A そのまま
//     k = 2: sqrtl(A)
//     k ∈ [3, 31]: std::pow(A, 1/k) を float で近似 → 整数化 → ±1 補正
//                  1/k は std::nextafter で「1/k 未満の最大の double」に丸めて
//                  std::pow が overshoot しないようにする
//     k ∈ [32, 63]: 答えは 1, 2, 3 のいずれか。 A vs 2^k と A vs 3^k で判定
//     k ≥ 64: A < 2^64 なので答えは 1
//
// per-query: ~10-20 cycle (std::pow が支配項、 k=2 は sqrtl で同様)
struct Solver {
 // overflow するなら 0 を返す
 static constexpr u64 sf_qpw(u64 a, int b, u64 r = 1) {
  for (; b;) {
   if (b & 1) { if (__builtin_umulll_overflow(a, r, &r)) return 0; }
   if (b >>= 1) { if (__builtin_umulll_overflow(a, a, &a)) return 0; }
  }
  return r;
 }

 // 注: constexpr にすると Apple clang で __builtin_umulll_overflow の compile-time
 // 評価が誤判定する (3^40 が overflow と判定される) のため runtime 評価に固定。
 template<size_t N>
 static auto make_ary(auto&& op) {
  std::array<decltype(op(0)), N> res{};
  for (size_t i = 0; i < N; ++i) res[i] = op(i);
  return res;
 }

 static u64 kth_root(u64 A, int k) {
  // 1/(i+3) で「1/k 未満の double」になるようにする (i=0,1,...,28 で k=3..31)
  static const auto iv = make_ary<29>([&](int i) { return std::nextafter(1.0 / (i + 3), 0.0); });
  // 3^k - 1 (k=32..63)。 A > pw3[k-32] なら答えが 4 以上 (起こらない)、
  // でも実装上は ans が 3 か 4 か (実は ≤ 3) を判定するため使用。
  static const auto pw3 = make_ary<32>([&](int i) { return sf_qpw(3, i + 32) - 1; });
  if (A == 0 || k == 1) return A;
  if (k < 32) {
   if (k == 2) {
    // ARM では long double = 64-bit double で sqrtl(2^64-1) が 2^32 に丸まる
    // (x64 の 80-bit long double では 4294967295.999... → 4294967295 で正しい)
    // 環境非依存にするため u128 で ±1 補正
    u64 r = u64(sqrtl((long double) A));
    while (r > 0 && u128(r) * r > A) --r;
    while (u128(r + 1) * (r + 1) <= A) ++r;
    return r;
   }
   u64 r = (i64) std::pow((double) A, iv[k - 3]);
   // (r+1)^k ≤ A なら r+1 が正解、 そうでなければ r。
   //   sf_qpw(r+1, k) - 1 < A   ↔   sf_qpw(r+1, k) ≤ A  (整数なので)
   return r + ((sf_qpw(r + 1, k) - 1) < A);
  }
  if (k > 63) return 1;
  // k ∈ [32, 63]: 答えは 1, 2, 3 のいずれか
  // ans ≥ 2 ⇔ A ≥ 2^k、 ans ≥ 3 ⇔ A > 3^k - 1 (= A ≥ 3^k)
  return 1 + (A >= (u64(1) << k)) + (A > pw3[k - 32]);
 }

 static std::vector<u64> run(const std::vector<std::pair<u64, int>>& queries) {
  std::vector<u64> ans;
  ans.reserve(queries.size());
  for (auto& [A, k] : queries) ans.push_back(kth_root(A, k));
  return ans;
 }
};
