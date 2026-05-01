#pragma once
#include "_common.hpp"
// =============================================================================
// Source: yosupo "Matrix Product (mod 2)" 提出 246187 を抽出。
//   https://judge.yosupo.jp/submission/246187
//   author: toxicpie (2024-10-16)
//   "I might have found the most cursed optimization technique in all of
//    competitive programming."
//   方式: 実行時に AVX2 アセンブリを JIT 生成して行列積を計算する。
//   各 12-row × 256-col ブロックごとに、A の bit パターンを基に VPXOR 命令の
//   レジスタ指定を書き換えた kernel を生成 (ビット = 1 の行だけ XOR される)。
//   JIT した kernel が B の列ブロックを舐めて C を書き換える。
// 抽出方針:
//   - I/O 部分 (read_input/write_output) は base.cpp 経由で string 入出力に
//   - main() は Conv 風 wrapper に
//   - x86_64 専用なので他環境では実行時に abort
//   - 全体を namespace yosupo_246187 で囲む
// ライセンス: 不明 (yosupo judge の慣習に従い参照、改造は最小限)
// =============================================================================

#ifndef __x86_64__
// 非 x86_64 環境では JIT/inline asm が使えないので runtime stub。
struct Mul {
 static vector<string> run(int, int, int, const vector<string>&, const vector<string>&) {
  fprintf(stderr, "yosupo_246187 (JIT): x86_64 required\n");
  std::abort();
 }
};
#else

#pragma GCC optimize("O3,unroll-loops,rename-registers")
#include <immintrin.h>
#include <sys/mman.h>
#include <unistd.h>
#include <cstring>

namespace yosupo_246187 {
using ::u32;
using ::u64;
using ::u8;

constexpr size_t N = 4096;
constexpr size_t LEN = N / 256;
constexpr size_t BLOCK_AR = 12;
constexpr size_t BLOCK_AC = 256;

[[gnu::aligned(64)]] inline __m256i A[N + BLOCK_AR][LEN + 1];
[[gnu::aligned(64)]] inline __m256i B[LEN][N + 74];
[[gnu::aligned(64)]] inline __m256i C[N + BLOCK_AR][LEN + 1];

typedef void (*code_ptr)(void*, void*);

constexpr size_t KERNEL_1_LEN = 128;
constexpr size_t KERNEL_2_LEN = 64;
constexpr size_t KERNEL_3_LEN = 93;

[[gnu::aligned(64)]] inline __m256i kernel_2_cache[KERNEL_2_LEN / 32];

[[gnu::naked]] inline void kernel_1() {
 __asm__(R"(
.intel_syntax noprefix
    vmovdqa ymm4,   YMMWORD PTR [rdi + 0x0000]
    vmovdqa ymm5,   YMMWORD PTR [rdi + 0x0220]
    vmovdqa ymm6,   YMMWORD PTR [rdi + 0x0440]
    vmovdqa ymm7,   YMMWORD PTR [rdi + 0x0660]
    vmovdqa ymm8,   YMMWORD PTR [rdi + 0x0880]
    vmovdqa ymm9,   YMMWORD PTR [rdi + 0x0aa0]
    vmovdqa ymm10,  YMMWORD PTR [rdi + 0x0cc0]
    vmovdqa ymm11,  YMMWORD PTR [rdi + 0x0ee0]
    vmovdqa ymm12,  YMMWORD PTR [rdi + 0x1100]
    vmovdqa ymm13,  YMMWORD PTR [rdi + 0x1320]
    vmovdqa ymm14,  YMMWORD PTR [rdi + 0x1540]
    vmovdqa ymm15,  YMMWORD PTR [rdi + 0x1760]
    mov     rdx,    0x20
    mov     rcx,    0x40
    vpxor   ymm0,   ymm0,   ymm0
    .nops 18
.att_syntax
    )");
}

[[gnu::naked]] inline void kernel_2() {
 __asm__(R"(
.intel_syntax noprefix
    vmovdqa ymm1,   YMMWORD PTR [rsi]
    vmovdqa ymm2,   YMMWORD PTR [rsi + rdx]
    vpxor   ymm3,   ymm1,   ymm2
    add     rsi,    rcx
    vpxor   ymm4,   ymm4,   ymm0
    vpxor   ymm5,   ymm5,   ymm0
    vpxor   ymm6,   ymm6,   ymm0
    vpxor   ymm7,   ymm7,   ymm0
    vpxor   ymm8,   ymm8,   ymm0
    vpxor   ymm9,   ymm9,   ymm0
    vpxor   ymm10,  ymm10,  ymm0
    vpxor   ymm11,  ymm11,  ymm0
    vpxor   ymm12,  ymm12,  ymm0
    vpxor   ymm13,  ymm13,  ymm0
    vpxor   ymm14,  ymm14,  ymm0
    vpxor   ymm15,  ymm15,  ymm0
.att_syntax
    )");
}

[[gnu::naked]] inline void kernel_3() {
 __asm__(R"(
.intel_syntax noprefix
    vmovdqa YMMWORD PTR [rdi + 0x0000], ymm4
    vmovdqa YMMWORD PTR [rdi + 0x0220], ymm5
    vmovdqa YMMWORD PTR [rdi + 0x0440], ymm6
    vmovdqa YMMWORD PTR [rdi + 0x0660], ymm7
    vmovdqa YMMWORD PTR [rdi + 0x0880], ymm8
    vmovdqa YMMWORD PTR [rdi + 0x0aa0], ymm9
    vmovdqa YMMWORD PTR [rdi + 0x0cc0], ymm10
    vmovdqa YMMWORD PTR [rdi + 0x0ee0], ymm11
    vmovdqa YMMWORD PTR [rdi + 0x1100], ymm12
    vmovdqa YMMWORD PTR [rdi + 0x1320], ymm13
    vmovdqa YMMWORD PTR [rdi + 0x1540], ymm14
    vmovdqa YMMWORD PTR [rdi + 0x1760], ymm15
    ret
.att_syntax
    )");
}

inline code_ptr jit_init() {
 code_ptr code = (code_ptr) mmap(NULL, 0x80000, PROT_READ | PROT_WRITE | PROT_EXEC, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
 u8* write_ptr = (u8*) code;
 std::memcpy((void*) write_ptr, (void*) kernel_1, KERNEL_1_LEN);
 std::memcpy(kernel_2_cache, (void*) kernel_2, KERNEL_2_LEN);
 std::memcpy((void*) (write_ptr + KERNEL_1_LEN + KERNEL_2_LEN * BLOCK_AC / 2), (void*) kernel_3, KERNEL_3_LEN);
 return code;
}

[[gnu::always_inline]] inline __m256i _mm256_expandmask1_2x16(uint16_t mask) {
 __m256i result = _mm256_set1_epi16(mask);
 __m256i shifts = _mm256_setr_epi32(0, 2, 4, 6, 8, 10, 12, 14);
 return _mm256_srlv_epi32(result, shifts);
}

[[gnu::always_inline]] inline __m256i _mm256_expandmask2_2x16(uint16_t mask) {
 __m256i result = _mm256_set1_epi16(mask);
 __m256i shifts = _mm256_setr_epi32(24, 22, 20, 18, 16, 14, 12, 10);
 __m256i select = _mm256_set1_epi32(0x03000000);
 return _mm256_sllv_epi32(result, shifts) & select;
}

inline void jit_compile_kernel(code_ptr code, size_t r, size_t c) {
 [[gnu::aligned(64)]] u32 masks[BLOCK_AC / 2] = {};
 for (size_t i = 0; i < BLOCK_AR; i++) {
  __m256i bits = _mm256_set1_epi32(0x03 << (i * 2 + 8));
  for (size_t j = 0; j < BLOCK_AC; j += 16) {
   uint16_t num = ((uint16_t*) A[r + i])[(c + j) / 16];
   ((__m256i*) masks)[j / 16] |= bits & _mm256_slli_epi32(_mm256_expandmask1_2x16(num), i * 2 + 8);
  }
 }
 for (size_t j = 0; j < BLOCK_AC / 2; j++) {
  __m256i* block_ptr = (__m256i*) ((u8*) code + KERNEL_1_LEN + KERNEL_2_LEN * j);
  for (size_t k = 0; k < 2; k++) {
   __m256i bits = _mm256_expandmask2_2x16(masks[j] >> (k * 16));
   __m256i vpxor = kernel_2_cache[k];
   _mm256_store_si256(&block_ptr[k], vpxor | bits);
  }
 }
}

} // namespace yosupo_246187

struct Mul {
 static vector<string> run(int n, int m, int k, const vector<string>& a, const vector<string>& b) {
  using namespace yosupo_246187;
  // 入力サイズが N=4096 を超えないことを前提。未使用領域は 0 で埋めるため毎回クリア。
  std::memset(A, 0, sizeof(A));
  std::memset(B, 0, sizeof(B));
  std::memset(C, 0, sizeof(C));
  // A: 行ごとに M ビットを 32-bit lane に詰める (((u32*)A[i])[j/32] の j%32 ビット目)。
  for (int i = 0; i < n; ++i) {
   const char* row = a[i].data();
   u32* p = (u32*) A[i];
   for (int j = 0; j < m; ++j)
    if (row[j] == '1') p[j / 32] |= u32(1) << (j % 32);
  }
  // B: 256-col ブロック単位の特殊レイアウト ((u32*)B[j/256])[i*8 + (j%256)/32] の j%32 ビット目。
  for (int i = 0; i < m; ++i) {
   const char* row = b[i].data();
   for (int j = 0; j < k; ++j)
    if (row[j] == '1')
     ((u32*) B[j / 256])[i * 8 + (j % 256) / 32] |= u32(1) << (j % 32);
  }
  static code_ptr code = nullptr;
  if (!code) code = jit_init();
  for (size_t c = 0; c < N; c += BLOCK_AC) {
   for (size_t r = 0; r < N; r += BLOCK_AR) {
    jit_compile_kernel(code, r, c);
    for (size_t i = 0; i < LEN; ++i) code(&C[r][i], &B[i][c]);
   }
  }
  vector<string> res(n, string(k, '0'));
  for (int i = 0; i < n; ++i) {
   const u32* src32 = (u32*) C[i];
   for (int j = 0; j < k; ++j)
    if ((src32[j / 32] >> (j % 32)) & 1) res[i][j] = '1';
  }
  return res;
 }
};

#endif // __x86_64__
