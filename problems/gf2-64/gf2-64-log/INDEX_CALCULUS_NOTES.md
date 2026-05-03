# gf2-64-log への Maspy 流 index calculus 適用 — 検討メモ

参考: https://maspypy.com/o1-mod-inv-mod-pow および
      https://judge.yosupo.jp/submission/354636 (purplesyringa)

## 構造の対応

| 概念 | Z/p^* (Maspy 元) | F_{2^64}^* (我々のターゲット) |
|---|---|---|
| 元の表現 | 整数 ≤ p ≈ 10^9 | F_2[X]/P(X) の多項式 (deg < 64) |
| Euclidean | 整数 GCD | F_2[X] の多項式 GCD |
| 「小さい」 | |a|, b ≤ √p ≈ 31623 | deg(a), deg(b) ≤ 32 |
| 因子基底 | 小素数 (≤ 31624 で ~3400 個) | 小次数既約多項式 (deg ≤ D で 2^D/D 個) |
| 因子基底 log 計算 | 1 回 BSGS / 各 = ~3ms × 3400 ≈ 10s | 同じく 1 回 BSGS / 各 |
| smoothness 確率 | ~30% | ~5-25% (deg ≤ 12-16 の factor base で) |
| Farey table | 1M bucket、bucket → (num/den) 直引き | **多項式版 Farey は実装煩雑** |
| Per query | rational approx + factor + lookup ≈ O(1) | 同等の実装が可能だが coding cost 高 |

## 主な実装段階

1. **F_2[X] の Euclidean / continued fraction** (~50 行)
   - x ∈ F_{2^64}^* に対し x · b(X) ≡ a(X) mod P(X) で deg(a), deg(b) ≤ 32
   - 標準的な拡張ユークリッドアルゴリズム

2. **F_2[X] の多項式因数分解** (~200 行)
   - Distinct-degree factorization: gcd(f, x^{2^k} - x) で degree-k 既約因子を抽出
   - Equal-degree factorization: Cantor-Zassenhaus (確率的)

3. **因子基底生成** (~100 行)
   - Rabin の既約性判定で deg ≤ D の既約多項式を列挙
   - F_2[X] 上の "small primes" 相当

4. **因子基底 log の事前計算** (~100 行)
   - 各既約多項式 q について log_g(q) mod (2^64-1) を BSGS で求める
   - 4 つの素数 (65535, 641, 65537, 6700417) ごとに CRT で分けて BSGS
   - **これが最大の init コスト** (~factor_base_size × ~30K cycle)

5. **Per-query 処理** (~50 行)
   - Euclidean → (a, b)
   - a, b を因子基底に factorize
   - smooth であれば log(x) = Σ log(factor_a) - Σ log(factor_b) mod (2^64-1)
   - smooth でなければ fallback (= 通常の BSGS)

合計: **~500 行**、テスト含めると ~700-1000 行。

## 性能試算

### 因子基底サイズ別の trade-off

| deg ≤ D | 因子基底サイズ | smoothness | init | per query (smooth) | per query (期待値) |
|---|---|---|---|---|---|
| 8 | ~30 | ~1% | 0.4 ms | 1K cycle | 30K cycle |
| 12 | ~335 | ~5% | 10 ms | 1K cycle | 28K cycle |
| 16 | ~4080 | ~25% | 120 ms | 1K cycle | 22K cycle |
| 20 | ~50000 | ~80% | 1500 ms | 1K cycle | 7K cycle |
| 32 | ~134M | 100% | (不可) | — | — |

**T による損益分岐**:

```
total_time(T, D) = init(D) + T × per_query(D)

T = 10K (我々のテスト):
   pohlig_v5: 0 + 10K × 13K cycle ≈ 130 ms
   indexcalc(deg ≤ 16): 120 + 10K × 22K = 120 + 220 = 340 ms (悪化)
   indexcalc(deg ≤ 20): 1500 + 10K × 7K = 1500 + 70 = 1570 ms (大悪化)

T = 100K:
   pohlig_v5: 1300 ms
   indexcalc(deg ≤ 16): 120 + 100K × 22K = 120 + 2200 = 2320 ms (悪化)
   indexcalc(deg ≤ 20): 1500 + 100K × 7K = 1500 + 700 = 2200 ms (悪化)

T = 1M:
   pohlig_v5: 13000 ms
   indexcalc(deg ≤ 16): 120 + 22000 = 22120 ms (悪化)
   indexcalc(deg ≤ 20): 1500 + 7000 = 8500 ms (改善 35%)

T = 10M:
   pohlig_v5: 130000 ms
   indexcalc(deg ≤ 20): 1500 + 70000 = 71500 ms (改善 45%)
   indexcalc(deg ≤ 24): もっと良い (smoothness ~95%)
```

つまり **T ≥ 1M 程度の超大量クエリでないと採算が合わない**。

## なぜ Z/p では効くが F_{2^64} では効きにくいのか

1. **smoothness 確率の差**:
   - Z 上 deg/log は連続的、Dickman ρ で smoothness 高い
   - F_2[X] 上の多項式は離散的、smoothness が指数的に低い

2. **Farey table の効率**:
   - Z/p の Farey table は 1M buckets の整数 indexing で per-query 1 lookup
   - F_2[X] は連分数展開を毎回 GCD で計算する必要、cache 化困難

3. **因子基底の基底単位**:
   - Z/p では小素数 < 32K で 3400 個の factor base 十分
   - F_2[X] では deg ≤ 16 でも 4080 個必要、smoothness ratio 低いので実質もっと欲しい

## 結論

- **本格移植は研究プロジェクト** (実装 ~1 週間)
- **我々のテスト規模 (T=10K-1M) では pohlig_v5 が事実上最適**
- T が極端に大きい (≥ 10M) なら Maspy 流が勝ちうる、その場合のみ実装価値あり

将来やるなら:
1. T ≥ 10M のテストケースを問題に追加
2. callgrind 用に T 100K + amortized init 込みで測定
3. 因子基底 log を Python script で offline 計算して header に embed (init 0)

## 関連実装

参考用に Z/p 版の Maspy 流を `yosupo-discrete-logarithm-fixed-mod/algos/index_calculus.hpp`
として移植済み。アルゴリズムの全体像はそちらで確認可能。
