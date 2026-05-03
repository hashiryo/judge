#pragma once
#ifdef USE_SIMDE
#include <simde/x86/avx2.h>
#endif
#include <bits/stdc++.h>
using namespace std;
using u8 = unsigned char;
using u32 = unsigned;
using i64 = long long;
using u64 = unsigned long long;
using u128 = __uint128_t;

// yosupo Stern-Brocot Tree (https://judge.yosupo.jp/problem/stern_brocot_tree)
//
// 入力フォーマット:
//   T 行のクエリ。各クエリは 5 種類のいずれか:
//   ENCODE_PATH a b
//     有理数 a/b を Stern-Brocot 木の根からの L/R パスとして出力。
//     path = "L k1 R k2 ..." の形式 (連続する同方向は連長圧縮)。
//   DECODE_PATH k1 k2 ... kN
//     L/R 連の長さ列から rational a/b を復元して出力。
//   LCA a1 b1 a2 b2
//     2 つの rational a1/b1, a2/b2 の LCA in SB tree、答えの rational を出力。
//   ANCESTOR depth a b
//     a/b の祖先で深さ depth のもの (= 一定数親に上がった先) を出力。 depth が
//     深すぎる場合は -1。
//   RANGE a b
//     a/b の Stern-Brocot interval (lx ly rx ry) を出力 (= mediant が a/b になる
//     左右の Stern-Brocot 隣接 fraction)。
//
// I/O が複雑なので harness は struct Solver::run() で全部処理させる。
