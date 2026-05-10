#pragma once
#include "_common.hpp"
// yosupo submission #351055 (https://judge.yosupo.jp/submission/351055) の移植。
// 当時 yosupo 3 番目の高速提出 (unique なアプローチ)。
//
// アプローチ: **Schönhage-Strassen radix-3** ですが #369907 と異なる実装。
//
// 核となるデータ構造 Circ3 = (v0, v1) は環 R[ω]/(ω² + ω + 1) を表現 (= F_2 上で
// ω² = ω + 1, 1 の原始 3 乗根)。 ω 倍を `cy1`, ω² 倍を `cy2` という rotate 演算
// で実装、 GF(2^64) 乗算なしで済む。
//
// 構造:
//   - Premult3 / Postmult3: 環の twist 変換
//   - FFT3: radix-3 FFT in Circ3
//   - CRT3: 2 つの並列計算 (Main3 を 2 回再帰呼び) を CRT で合成
//   - Main3: l × m 分割で再帰、 各 l ブロックを Main3(...) で処理
//
// #369907 (schoenhage_radix3) との違い:
//   - 環の選択: schoenhage_radix3 は次元 doubling (size 2d)、 schoenhage_strassen3 は
//     Circ3 (2 component) で 2 並列 (c1, c2) を CRT 合成
//   - butterfly の twiddle はどちらも cyclic shift (= rotate + XOR)
//
// 必要拡張: PCLMUL + SSE4.1。

#pragma GCC optimize("O3,unroll-loops")

#include <algorithm>
#include <vector>
namespace conv_f2_64_schoenhage_strassen3 {
// =============================================================================
// F2XMod1000000000000001B: GF(2^64) = F_2[x] / (x^64 + x^4 + x^3 + x + 1)
// (元コードの x*x → 同一引数 clmul の最適化付き mul)
// =============================================================================
struct F2XMod1000000000000001B {
 using uu64= unsigned long long;
 uu64 x;
 PCLMUL static uu64 mul(uu64 a, uu64 b) {
  __m128i ab= _mm_set_epi64x(a, b);
  __m128i xy= _mm_clmulepi64_si128(ab, ab, 1);
  uu64 X= _mm_extract_epi64(xy, 0), Y= _mm_extract_epi64(xy, 1);
  X^= Y ^ Y << 1 ^ Y << 3 ^ Y << 4;
  Y= Y >> 63 ^ Y >> 61 ^ Y >> 60;
  X^= Y ^ Y << 1 ^ Y << 3 ^ Y << 4;
  return X;
 }
 using Self= F2XMod1000000000000001B;
 F2XMod1000000000000001B(uu64 v= 0): x(v) {}
 Self operator+(const Self& r) const { return Self(x ^ r.x); }
 Self operator-(const Self& r) const { return Self(x ^ r.x); }
 Self operator*(const Self& r) const { return Self(mul(x, r.x)); }
 Self& operator+=(const Self& r) {
  x^= r.x;
  return *this;
 }
 Self& operator-=(const Self& r) {
  x^= r.x;
  return *this;
 }
 Self& operator*=(const Self& r) {
  x= mul(x, r.x);
  return *this;
 }
 Self operator-() const { return Self(x); }
};
using RingElem= F2XMod1000000000000001B;
// =============================================================================
// IsThisSchonhageStrassen3 本体 (元コードそのまま、 ring_mul_operation_count は除去)
// =============================================================================
struct IsThisSchonhageStrassen3 {
 struct Circ3 {
  RingElem v0, v1;
  Circ3 cy1() const { return {-v1, v0 - v1}; }
  Circ3 cy2() const { return {v1 - v0, -v0}; }
  Circ3& operator+=(const Circ3& r) {
   v0+= r.v0;
   v1+= r.v1;
   return *this;
  }
  Circ3& operator-=(const Circ3& r) {
   v0-= r.v0;
   v1-= r.v1;
   return *this;
  }
  Circ3 operator+(const Circ3& r) const { return {v0 + r.v0, v1 + r.v1}; }
  Circ3 operator-(const Circ3& r) const { return {v0 - r.v0, v1 - r.v1}; }
  Circ3 operator*(const Circ3& r) const {
   RingElem a= v0 * r.v0, b= v1 * r.v1;
   return {a - b, (v0 + v1) * (r.v0 + r.v1) - a - b - b};
  }
  Circ3 rev() const { return {v0 - v1, -v1}; }
 };
 RingElem zero;
 Circ3 zero3;
 IsThisSchonhageStrassen3(RingElem _zero): zero(_zero), zero3({zero, zero}) {}
 void mulAdd3(std::vector<Circ3>::iterator a, std::vector<Circ3>::iterator b, int m, int sh, int d) {
  if(sh == 0) {
   for(int i= 0; i < m - d; i++) b[i + d]+= a[i];
   for(int i= m - d; i < m; i++) b[i - m + d]+= a[i].cy1();
  } else if(sh == 1) {
   for(int i= 0; i < m - d; i++) b[i + d]+= a[i].cy1();
   for(int i= m - d; i < m; i++) b[i - m + d]+= a[i].cy2();
  } else {
   for(int i= 0; i < m - d; i++) b[i + d]+= a[i].cy2();
   for(int i= m - d; i < m; i++) b[i - m + d]+= a[i];
  }
 }
 void mulAdd3(std::vector<Circ3>::iterator a, std::vector<Circ3>::iterator b, int m, int i) { mulAdd3(a, b, m, i / m % 3, i % m); }
 void Premult3(std::vector<Circ3>& a, int m, int l) {
  int gap= m / l;
  std::vector<Circ3> b(m * l, zero3);
  for(int i= 0; i < l; i++) mulAdd3(a.begin() + i * m, b.begin() + i * m, m, gap * i);
  std::swap(a, b);
 }
 void Postmult3(std::vector<Circ3>& a, int m, int l) {
  int gap= m / l;
  std::vector<Circ3> b(m * l, zero3);
  for(int i= 0; i < l; i++) mulAdd3(a.begin() + i * m, b.begin() + i * m, m, m * 3 - gap * i);
  std::swap(a, b);
 }
 void FFT3(std::vector<Circ3>& a, int m, int l) {
  int gap= m / l, msdig= l / 3;
  for(int z= l / 3; z; z/= 3) {
   std::vector<Circ3> b(a.size(), zero3);
   for(int p= 0; p < msdig; p+= z)
    for(int q= 0; q < z; q++)
     for(int t= 0; t < 3; t++)
      for(int u= 0; u < 3; u++) {
       mulAdd3(a.begin() + ((p * 3 + u * z + q) * m), b.begin() + (p + q + msdig * t) * m, m, (l * t * u + p * u * 3) * gap);
      }
   std::swap(a, b);
  }
 }
 std::pair<Circ3, Circ3> CRT3(Circ3 c1, Circ3 c2) {
  Circ3 d= c2 - c1;
  d= d + d.cy1() + d.cy1();
  return {c1 + c1 + c1 - d.cy1(), d};
 }
 void ElemMult3(std::vector<Circ3>::iterator a, std::vector<Circ3>::iterator b, std::vector<Circ3>::iterator c) {
  Circ3 u= a[0] * b[0], v= a[1] * b[1], w= a[2] * b[2];
  Circ3 uv= (a[0] + a[1]) * (b[0] + b[1]), vw= (a[1] + a[2]) * (b[1] + b[2]);
  c[0]= u + (vw - v - w).cy1();
  c[1]= uv - u - v + w.cy1();
  c[2]= (a[0] + a[1] + a[2]) * (b[0] + b[1] + b[2]) - uv - vw + v + v;
 }
 void Main3(std::vector<Circ3> a, std::vector<Circ3> b, std::vector<Circ3>::iterator dest, int logn) {
  if(logn == 1) {
   ElemMult3(a.begin(), b.begin(), dest);
  } else {
   int logl= logn / 2, logm= logn - logl;
   int l= 1;
   for(int i= 0; i < logl; i++) l*= 3;
   int m= 1;
   for(int i= 0; i < logm; i++) m*= 3;
   int n= l * m;
   std::vector<Circ3> a2(n, zero3), b2(n, zero3), c1(n, zero3), c2(n, zero3);
   for(int i= 0; i < n; i++) a2[i]= a[i].rev();
   for(int i= 0; i < n; i++) b2[i]= b[i].rev();
   Premult3(a, m, l);
   Premult3(b, m, l);
   Postmult3(a2, m, l);
   Postmult3(b2, m, l);
   FFT3(a, m, l);
   FFT3(b, m, l);
   FFT3(a2, m, l);
   FFT3(b2, m, l);
   for(int i= 0; i < l; i++) {
    Main3(std::vector<Circ3>(a.begin() + i * m, a.begin() + (i * m + m)), std::vector<Circ3>(b.begin() + i * m, b.begin() + (i * m + m)), c1.begin() + i * m, logm);
    Main3(std::vector<Circ3>(a2.begin() + i * m, a2.begin() + (i * m + m)), std::vector<Circ3>(b2.begin() + i * m, b2.begin() + (i * m + m)), c2.begin() + i * m, logm);
   }
   FFT3(c1, m, l);
   for(int i= 1; i < l - i; i++)
    for(int j= 0; j < m; j++) std::swap(c1[i * m + j], c1[(l - i) * m + j]);
   FFT3(c2, m, l);
   for(int i= 1; i < l - i; i++)
    for(int j= 0; j < m; j++) std::swap(c2[i * m + j], c2[(l - i) * m + j]);
   Postmult3(c1, m, l);
   Premult3(c2, m, l);
   for(auto& x: c2) x= x.rev();
   std::vector<Circ3> c(a.size(), zero3);
   for(int i= 0; i < (l - 1) * m; i++) {
    auto [u, v]= CRT3(c1[i], c2[i]);
    c[i]+= u;
    c[i + m]+= v;
   }
   for(int i= (l - 1) * m; i < l * m; i++) {
    auto [u, v]= CRT3(c1[i], c2[i]);
    c[i]+= u;
    c[i + m - n]+= v.cy1();
   }
   std::copy(c.begin(), c.end(), dest);
  }
 }
};
inline std::vector<RingElem> IsThisSchonhageStrassenConvolution3(std::vector<RingElem> a, std::vector<RingElem> b, RingElem zero) {
 int target_n= (int)std::max(a.size(), b.size());
 int logn= 0;
 int n= 1;
 while(n < target_n) {
  logn++;
  n*= 3;
 }
 if(n == 1) {
  return {a[0] * b[0]};
 }
 auto ds= IsThisSchonhageStrassen3(zero);
 IsThisSchonhageStrassen3::Circ3 zero3= {zero, zero};
 std::vector<IsThisSchonhageStrassen3::Circ3> a2(n, zero3), b2(n, zero3), c2(n, zero3);
 for(int i= 0; i < (int)a.size(); i++) a2[i].v0= a[i];
 for(int i= 0; i < (int)b.size(); i++) b2[i].v0= b[i];
 ds.Main3(a2, b2, c2.begin(), logn);
 int csz= a.size() + b.size() - 1;
 std::vector<RingElem> c(csz, {0});
 for(int i= 0; i < csz; i++) c[i]= i < n ? c2[i].v0 : c2[i - n].v1;
 return c;
}
}  // namespace conv_f2_64_schoenhage_strassen3
struct Solver {
 static std::vector<u64> run(int n, int m, const std::vector<u64>& a_in, const std::vector<u64>& b_in) {
  using namespace conv_f2_64_schoenhage_strassen3;
  std::vector<RingElem> a(n), b(m);
  for(int i= 0; i < n; ++i) a[i]= RingElem(a_in[i]);
  for(int i= 0; i < m; ++i) b[i]= RingElem(b_in[i]);
  auto c= IsThisSchonhageStrassenConvolution3(std::move(a), std::move(b), RingElem(0));
  std::vector<u64> out(c.size());
  for(size_t i= 0; i < c.size(); ++i) out[i]= c[i].x;
  return out;
 }
};
