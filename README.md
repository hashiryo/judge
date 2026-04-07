# judge

競プロの作業用ジャッジシステム。問題単位で複数の提出を比較できる。

## 使い方

`problems/` 以下に問題ディレクトリを作り、`.cpp` ファイルを置いて push する。

```
problems/
  yosupo-unionfind/
    sol_a.cpp
    sol_b.cpp
```

各 `.cpp` ファイルの先頭に問題 URL を記載:

```cpp
// judge: https://judge.yosupo.jp/problem/unionfind
```

自作テストケースを使う場合は `testcases/` ディレクトリに配置:

```
problems/my-problem/
  sol.cpp
  testcases/
    00.in  00.out
    01.in  01.out
```

push すると GHA が自動でテストを実行し、同じ問題の全提出を比較した結果を Step Summary に表示する。
