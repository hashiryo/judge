#pragma once
// =============================================================================
// 後方互換用の meta-header。 中身は分割されて mul.hpp / sq.hpp / pow.hpp /
// inv.hpp / sqrt.hpp に移動済み。
//
// 新しい algo は必要な分だけ個別 include することを推奨:
//   #include "../../_shared/mul.hpp"  ← mul だけ要る場合
//   #include "../../_shared/sq.hpp"   ← mul + sq 要る場合 (sq.hpp は mul.hpp を pull-in)
//   #include "../../_shared/pow.hpp"  ← pow まで要る場合
//   ...
//
// 既存の `#include "../../_shared/pclmul_core.hpp"` を直接 split-headers に変えても
// 動作は変わらない。 リファクタリングの一環で個別ファイルへ移行を推奨。
//
// 運用ポリシーは各 split file 冒頭のコメントを参照。
// =============================================================================
#include "mul.hpp"
#include "sq.hpp"
#include "pow.hpp"
#include "inv.hpp"
#include "sqrt.hpp"
