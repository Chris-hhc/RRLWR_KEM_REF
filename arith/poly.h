/* 
 * Copyright 2026 NXP
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef POLY_H
#define POLY_H

#include "fprime.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct{
    int32_t coeffs[RRLWR_N];
  } poly;

  static inline int32_t mod_q13_u(int64_t x) {
    return (int32_t)((uint64_t)x & RRLWR_Q_MASK);
  }

  static inline int32_t mod_q13_centered(int64_t x) {
    int32_t y = mod_q13_u(x);
    if(y > (RRLWR_Q >> 1)) {
      y -= RRLWR_Q;
    }
    return y;
  }

  static inline int32_t q13_reduce_u_i64(int64_t x) {
    return (int32_t)((uint64_t)x & RRLWR_Q_MASK);
  }

  static inline int32_t q13_center_i32(int32_t x) {
    int32_t y = x & RRLWR_Q_MASK;
    if(y > (RRLWR_Q >> 1)) {
      y -= RRLWR_Q;
    }
    return y;
  }

  static inline int32_t q13_small_secret_i32(int32_t x) {
    return q13_center_i32(x);
  }

  void poly_ntt32(poly *f, int32_t prime, int32_t primeinv, int32_t fp_zetas[RRLWR_N]);
  void poly_invntt32(poly *f, int32_t prime, int32_t primeinv, int32_t finalconst, int32_t fp_zetas[RRLWR_N]);
  void poly_basemul32(poly *r, const poly *f, const poly *g, int32_t prime, int32_t primeinv);
  void poly_mul_yplus2_q13(poly *r, const poly *a);
  void poly_mul_q13_schoolbook(poly *r, const poly *a, const poly *b);
  void poly_mul_q13_toom4x32_karatsuba(poly *r, const poly *a, const poly *b);
  void poly_mul_q13(poly *r, const poly *a, const poly *b);
  void poly_macc_q13(poly *acc, const poly *a, const poly *b);
  void poly_macc_q13_smallsecret(poly *acc, const poly *dense, const poly *small);
  void poly_add32(poly *r, const poly *f, const poly *g, int32_t prime);
  void poly_add(poly *r, poly *f, poly *g);
  void poly_sub32(poly *r, poly *f, poly *g, int32_t prime);
  void poly_sub(poly *r, poly *f, poly *g);
  void poly_reduce_pow2(poly *r, poly *f, int32_t d);
  void poly_conditional_final_reduce32(poly *r, int32_t prime);
  void poly_round_xtoy(poly *r, const poly *f, int32_t x, int32_t y);
  void poly_compress(poly *r, int32_t x);
  void poly_decompress(poly *r, int32_t x);

#ifdef __cplusplus
}
#endif

#endif
