#pragma once
#include "_common.hpp"
// 教科書的 Warshall-Floyd。3 重ループ、scalar。
// d[i*Vp + j] でアクセス。
struct WF {
 int V, Vp;
 i32* d;
 WF(int v, int vp, i32* dist): V(v), Vp(vp), d(dist) {}
 void run() {
  for (int k= 0; k < V; ++k) {
   for (int i= 0; i < V; ++i) {
    i32 dik= d[i * Vp + k];
    for (int j= 0; j < V; ++j) {
     i32 nd= dik + d[k * Vp + j];
     if (d[i * Vp + j] > nd) d[i * Vp + j]= nd;
    }
   }
  }
 }
};
