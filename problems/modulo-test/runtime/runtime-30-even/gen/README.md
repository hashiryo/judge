# runtime-30-even

入力で法 `mod` を与える、32 bit 動的 mod 用のベンチ問題。`mod` は **偶数** に制約される。

反復は

`s <- (a * s + b) * (c * s + d) (mod mod)`

で、`mod` は毎ケースごとに与えられる。

通常の Montgomery reduction は奇数 mod を要求するため、本問題では `mod = 2^s * mo` (mo は奇数) の形に分解して reduction する `montgomery_even*` 系統と、汎用 (Barrett / div2by1 / long_double / naive) を比較する。

## 入力

1 行に 7 整数を与える。

`N mod s a b c d`

- `0 <= N`
- `2^20 < mod < 2^30`
- `mod` は偶数
- `0 <= s, a, b, c, d < mod`

## 出力

`N` 回反復した後の `s` を出力する。

## ケース方針

- 小さい手書きケース
- 境界寄りの `mod`
- `2^20` 付近の最小側
- `2^30` 付近の最大側 / 2 のべき
- ランダムな偶数 mod

## 再生成

以下を問題ディレクトリで実行する。

```bash
./gen/generate.sh
```
