



- 現在最速は tastest_full_dif_radix4_v2_loctbl.hpp
- submask table を constexpr にしても工夫すれば loctbl 同等のスピードにならないかなあ.
  - sub_n は 一瞬で求められるでしょ.
  - ポインタで取得しなくていいのでは？
- 4基底なら reduction を遅延できる可能性があるかも？
  - (思いつき・未検証. 検証するなら fastest_simple_dif_radix4.hpp で考えるのが良いと思う. vmul 使っていないのでみやすいはず)
- vmul の適用の仕方が少しバラバラだから揃えてもいいかも？そんなことない？
