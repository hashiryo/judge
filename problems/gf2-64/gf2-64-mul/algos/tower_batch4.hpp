#pragma once
// Nimber 塔乗算を 4-way 並べて呼び、コンパイラの命令スケジューラに ILP を期待する
// バリアント。Nimber::mul の 11 個の F_{2^16} log/exp lookup を、4 入力分隣接させて
// 配置する (= 44 個の独立 lookup を OoO に晒す)。
//
// 実装上は Nimber::operator* を 4 つループで連続して呼ぶだけ。コンパイラが
// 全部 inline 展開してくれれば、各 step を batch4 で interleave する余地が
// 生まれる (実際には別 mul 同士の data dependency はゼロなので OoO に任せる
// だけでも十分なはずだが、register allocation のヒントとして手動 unroll する)。
#pragma GCC optimize("O3,unroll-loops")
#include "_common.hpp"
#include "../../_shared/basis_change.hpp"
#include "mylib/algebra/Nimber.hpp"

struct GF2_64Op {
 static vector<u64> run(const vector<u64>& as, const vector<u64>& bs) {
  using gf2_64_basis::poly_to_nim;
  using gf2_64_basis::nim_to_poly;
  Nimber::init();
  const size_t T = as.size();
  vector<u64> ans(T);
  size_t i = 0;
  for (; i + 4 <= T; i += 4) {
   Nimber A0(poly_to_nim(as[i  ])), B0(poly_to_nim(bs[i  ]));
   Nimber A1(poly_to_nim(as[i+1])), B1(poly_to_nim(bs[i+1]));
   Nimber A2(poly_to_nim(as[i+2])), B2(poly_to_nim(bs[i+2]));
   Nimber A3(poly_to_nim(as[i+3])), B3(poly_to_nim(bs[i+3]));
   // 4 個の独立 mul を並べる (互いに data dependency なし → OoO に並列化されるはず)
   Nimber C0 = A0 * B0;
   Nimber C1 = A1 * B1;
   Nimber C2 = A2 * B2;
   Nimber C3 = A3 * B3;
   ans[i  ] = nim_to_poly(C0.val());
   ans[i+1] = nim_to_poly(C1.val());
   ans[i+2] = nim_to_poly(C2.val());
   ans[i+3] = nim_to_poly(C3.val());
  }
  for (; i < T; ++i) {
   ans[i] = nim_to_poly((Nimber(poly_to_nim(as[i])) * Nimber(poly_to_nim(bs[i]))).val());
  }
  return ans;
 }
};
