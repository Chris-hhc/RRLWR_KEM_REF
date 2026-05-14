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

#include "poly.h"

#ifndef NDEBUG
#include <stdio.h>
#include <stdlib.h>
#endif

#define TOOM4_BLK (32)
#define TOOM4_BASE (16)
#define TOOM4_SMALL_PROD (31)
#define TOOM4_MID_PROD (63)
#define TOOM4_FULL_PROD (255)

static inline int32_t q13_reduce_u(int64_t x) {
  return (int32_t)((uint64_t)x & 8191u);
}

static inline int64_t q13_center_i64(int32_t x) {
  int32_t y = x & 8191;

  if(y > 4096) {
    y -= 8192;
  }

  return (int64_t)y;
}

static inline int64_t divexact_i64(int64_t x, int64_t d) {
#ifndef NDEBUG
  if(x % d != 0) {
    fprintf(stderr, "divexact_i64 failed: %lld %% %lld != 0\n",
            (long long)x, (long long)d);
    abort();
  }
#endif

  return x / d;
}

static void mul16_schoolbook_i64(int64_t out[TOOM4_SMALL_PROD],
                                 const int64_t a[TOOM4_BASE],
                                 const int64_t b[TOOM4_BASE]) {
  for(unsigned int i = 0; i < TOOM4_SMALL_PROD; i++) {
    out[i] = 0;
  }

  for(unsigned int i = 0; i < TOOM4_BASE; i++) {
    for(unsigned int j = 0; j < TOOM4_BASE; j++) {
      out[i + j] += a[i] * b[j];
    }
  }
}

static void mul32_karatsuba_i64(int64_t out[TOOM4_MID_PROD],
                                const int64_t a[TOOM4_BLK],
                                const int64_t b[TOOM4_BLK]) {
  int64_t z0[TOOM4_SMALL_PROD], z1[TOOM4_SMALL_PROD], z2[TOOM4_SMALL_PROD];
  int64_t sa[TOOM4_BASE], sb[TOOM4_BASE];

  for(unsigned int i = 0; i < TOOM4_MID_PROD; i++) {
    out[i] = 0;
  }

  mul16_schoolbook_i64(z0, a, b);
  mul16_schoolbook_i64(z2, a + TOOM4_BASE, b + TOOM4_BASE);

  for(unsigned int i = 0; i < TOOM4_BASE; i++) {
    sa[i] = a[i] + a[i + TOOM4_BASE];
    sb[i] = b[i] + b[i + TOOM4_BASE];
  }

  mul16_schoolbook_i64(z1, sa, sb);

  for(unsigned int i = 0; i < TOOM4_SMALL_PROD; i++) {
    z1[i] -= z0[i] + z2[i];
    out[i] += z0[i];
    out[i + TOOM4_BASE] += z1[i];
    out[i + 2 * TOOM4_BASE] += z2[i];
  }
}

void poly_ntt32(poly *f, int32_t prime, int32_t primeinv, int32_t fp_zetas[RRLWR_N]) {
  unsigned int len, start, j, k;
  int32_t t, fp_zeta;
  int32_t *fc = f->coeffs;

  k = 1;
  for(len = (RRLWR_N >> 1); len >= 1; len >>= 1) {
    for(start = 0; start < RRLWR_N; start = j + len) {
      fp_zeta = fp_zetas[k++];
      for(j = start; j < start + len; j++) {
        t = montgomery_mul32(fp_zeta, fc[j + len], prime, primeinv);
        fc[j] = conditional_reduce32(fc[j], prime); // Reduce back to [-p, p]
        fc[j + len] = fc[j] - t;
        fc[j] = fc[j] + t;
      }
    }
  }

  for(j = 0; j < RRLWR_N; j++) {
      fc[j] = conditional_reduce32(fc[j], prime); // Reduce back to [-p, p]
  }
}

void poly_invntt32(poly *f, int32_t prime, int32_t primeinv, int32_t finalconst, int32_t fp_zetas[RRLWR_N]) {
  unsigned int start, len, j, k;
  int32_t t, zeta;
  int32_t *fc = f->coeffs;

  k = RRLWR_N-1;
  for(len = 1; len <= (RRLWR_N >> 1); len <<= 1) {
    for(start = 0; start < RRLWR_N; start = j + len) {
      zeta = fp_zetas[k--];
      for(j = start; j < start + len; j++) {
        t = fc[j];
        fc[j] = t + fc[j + len];
        fc[j + len] = fc[j + len] - t;
        fc[j + len] = montgomery_mul32(zeta, fc[j + len], prime, primeinv);
        fc[j] = conditional_reduce32(fc[j], prime); // Reduce back to [-p, p]
      }
    }
  }

  for(j = 0; j < RRLWR_N; j++) {
    fc[j] = montgomery_mul32(fc[j], finalconst, prime, primeinv);
  }
}

void poly_basemul32(poly *r, const poly *f, const poly *g, int32_t prime, int32_t primeinv) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = montgomery_mul32(f->coeffs[i], g->coeffs[i], prime, primeinv);
  }
}

void poly_mul_yplus2_q13(poly *r, const poly *a) {
  r->coeffs[0] = mod_q13_u(-(int64_t)a->coeffs[RRLWR_N - 1] + 2LL * a->coeffs[0]);

  for(unsigned int i = 1; i < RRLWR_N; i++) {
    r->coeffs[i] = mod_q13_u((int64_t)a->coeffs[i - 1] + 2LL * a->coeffs[i]);
  }
}

void poly_mul_q13_schoolbook(poly *r, const poly *a, const poly *b) {
  int64_t tmp[RRLWR_N] = {0};

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    int64_t ai = a->coeffs[i];

    for(unsigned int j = 0; j < RRLWR_N; j++) {
      int64_t prod = ai * (int64_t)b->coeffs[j];
      unsigned int d = i + j;

      if(d < RRLWR_N) {
        tmp[d] += prod;
      } else {
        tmp[d - RRLWR_N] -= prod;
      }
    }
  }

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = mod_q13_u(tmp[i]);
  }
}

void poly_mul_q13_toom4x32_karatsuba(poly *r, const poly *a, const poly *b) {
  int64_t a0[TOOM4_BLK], a1[TOOM4_BLK], a2[TOOM4_BLK], a3[TOOM4_BLK];
  int64_t b0[TOOM4_BLK], b1[TOOM4_BLK], b2[TOOM4_BLK], b3[TOOM4_BLK];
  int64_t ae0[TOOM4_BLK], ae1[TOOM4_BLK], aem1[TOOM4_BLK];
  int64_t ae2[TOOM4_BLK], aem2[TOOM4_BLK], ae3[TOOM4_BLK], aeinf[TOOM4_BLK];
  int64_t be0[TOOM4_BLK], be1[TOOM4_BLK], bem1[TOOM4_BLK];
  int64_t be2[TOOM4_BLK], bem2[TOOM4_BLK], be3[TOOM4_BLK], beinf[TOOM4_BLK];
  int64_t v0[TOOM4_MID_PROD], v1[TOOM4_MID_PROD], vm1[TOOM4_MID_PROD];
  int64_t v2[TOOM4_MID_PROD], vm2[TOOM4_MID_PROD], v3[TOOM4_MID_PROD];
  int64_t vinf[TOOM4_MID_PROD];
  int64_t full[TOOM4_FULL_PROD];

  for(unsigned int i = 0; i < TOOM4_BLK; i++) {
    a0[i] = q13_center_i64(a->coeffs[i]);
    a1[i] = q13_center_i64(a->coeffs[TOOM4_BLK + i]);
    a2[i] = q13_center_i64(a->coeffs[2 * TOOM4_BLK + i]);
    a3[i] = q13_center_i64(a->coeffs[3 * TOOM4_BLK + i]);

    b0[i] = q13_center_i64(b->coeffs[i]);
    b1[i] = q13_center_i64(b->coeffs[TOOM4_BLK + i]);
    b2[i] = q13_center_i64(b->coeffs[2 * TOOM4_BLK + i]);
    b3[i] = q13_center_i64(b->coeffs[3 * TOOM4_BLK + i]);

    ae0[i] = a0[i];
    ae1[i] = a0[i] + a1[i] + a2[i] + a3[i];
    aem1[i] = a0[i] - a1[i] + a2[i] - a3[i];
    ae2[i] = a0[i] + 2 * a1[i] + 4 * a2[i] + 8 * a3[i];
    aem2[i] = a0[i] - 2 * a1[i] + 4 * a2[i] - 8 * a3[i];
    ae3[i] = a0[i] + 3 * a1[i] + 9 * a2[i] + 27 * a3[i];
    aeinf[i] = a3[i];

    be0[i] = b0[i];
    be1[i] = b0[i] + b1[i] + b2[i] + b3[i];
    bem1[i] = b0[i] - b1[i] + b2[i] - b3[i];
    be2[i] = b0[i] + 2 * b1[i] + 4 * b2[i] + 8 * b3[i];
    bem2[i] = b0[i] - 2 * b1[i] + 4 * b2[i] - 8 * b3[i];
    be3[i] = b0[i] + 3 * b1[i] + 9 * b2[i] + 27 * b3[i];
    beinf[i] = b3[i];
  }

  mul32_karatsuba_i64(v0, ae0, be0);
  mul32_karatsuba_i64(v1, ae1, be1);
  mul32_karatsuba_i64(vm1, aem1, bem1);
  mul32_karatsuba_i64(v2, ae2, be2);
  mul32_karatsuba_i64(vm2, aem2, bem2);
  mul32_karatsuba_i64(v3, ae3, be3);
  mul32_karatsuba_i64(vinf, aeinf, beinf);

  for(unsigned int i = 0; i < TOOM4_FULL_PROD; i++) {
    full[i] = 0;
  }

  for(unsigned int t = 0; t < TOOM4_MID_PROD; t++) {
    int64_t c0 = v0[t];
    int64_t c6 = vinf[t];
    int64_t c1 = divexact_i64(-40 * v0[t] + 120 * v1[t] - 60 * vm1[t]
                              - 30 * v2[t] + 6 * vm2[t] + 4 * v3[t]
                              - 1440 * vinf[t], 120);
    int64_t c2 = divexact_i64(-150 * v0[t] + 80 * v1[t] + 80 * vm1[t]
                              - 5 * v2[t] - 5 * vm2[t]
                              + 480 * vinf[t], 120);
    int64_t c3 = divexact_i64(50 * v0[t] - 70 * v1[t] - 5 * vm1[t]
                              + 35 * v2[t] - 5 * vm2[t] - 5 * v3[t]
                              + 1800 * vinf[t], 120);
    int64_t c4 = divexact_i64(30 * v0[t] - 20 * v1[t] - 20 * vm1[t]
                              + 5 * v2[t] + 5 * vm2[t]
                              - 600 * vinf[t], 120);
    int64_t c5 = divexact_i64(-10 * v0[t] + 10 * v1[t] + 5 * vm1[t]
                              - 5 * v2[t] - vm2[t] + v3[t]
                              - 360 * vinf[t], 120);

    full[t] += c0;
    full[TOOM4_BLK + t] += c1;
    full[2 * TOOM4_BLK + t] += c2;
    full[3 * TOOM4_BLK + t] += c3;
    full[4 * TOOM4_BLK + t] += c4;
    full[5 * TOOM4_BLK + t] += c5;
    full[6 * TOOM4_BLK + t] += c6;
  }

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    int64_t x = full[i];

    if(i + RRLWR_N < TOOM4_FULL_PROD) {
      x -= full[i + RRLWR_N];
    }

    r->coeffs[i] = q13_reduce_u(x);
  }
}

void poly_mul_q13(poly *r, const poly *a, const poly *b) {
  poly_mul_q13_toom4x32_karatsuba(r, a, b);
}

void poly_macc_q13(poly *acc, const poly *a, const poly *b) {
  poly t;

  poly_mul_q13(&t, a, b);

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    acc->coeffs[i] = mod_q13_u((int64_t)acc->coeffs[i] + t.coeffs[i]);
  }
}

void poly_macc_q13_smallsecret(poly *acc, const poly *dense, const poly *small) {
  int64_t tmp[RRLWR_N];
  int32_t dense_c[RRLWR_N];

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    tmp[i] = q13_center_i32(acc->coeffs[i]);
    dense_c[i] = q13_center_i32(dense->coeffs[i]);
  }

  /*
   * This ref fast path branches on small secret coefficients; replace with
   * constant-time masks if required.
   */
  for(unsigned int c = 0; c < RRLWR_N; c++) {
    int32_t sc = q13_small_secret_i32(small->coeffs[c]);

#ifndef NDEBUG
    if(!(sc == -2 || sc == -1 || sc == 0 || sc == 1 || sc == 2)) {
      fprintf(stderr, "poly_macc_q13_smallsecret: non-small coeff %d\n", sc);
      abort();
    }
#endif

    if(sc == 0) {
      continue;
    }

    if(sc == 1) {
      for(unsigned int t = 0; t < RRLWR_N - c; t++) {
        tmp[t + c] += dense_c[t];
      }
      for(unsigned int t = RRLWR_N - c; t < RRLWR_N; t++) {
        tmp[t + c - RRLWR_N] -= dense_c[t];
      }
    } else if(sc == -1) {
      for(unsigned int t = 0; t < RRLWR_N - c; t++) {
        tmp[t + c] -= dense_c[t];
      }
      for(unsigned int t = RRLWR_N - c; t < RRLWR_N; t++) {
        tmp[t + c - RRLWR_N] += dense_c[t];
      }
    } else if(sc == 2) {
      for(unsigned int t = 0; t < RRLWR_N - c; t++) {
        tmp[t + c] += 2LL * dense_c[t];
      }
      for(unsigned int t = RRLWR_N - c; t < RRLWR_N; t++) {
        tmp[t + c - RRLWR_N] -= 2LL * dense_c[t];
      }
    } else {
      for(unsigned int t = 0; t < RRLWR_N - c; t++) {
        tmp[t + c] -= 2LL * dense_c[t];
      }
      for(unsigned int t = RRLWR_N - c; t < RRLWR_N; t++) {
        tmp[t + c - RRLWR_N] += 2LL * dense_c[t];
      }
    }
  }

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    acc->coeffs[i] = q13_reduce_u_i64(tmp[i]);
  }
}

void poly_add32(poly *r, const poly *f, const poly *g, int32_t prime) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = add32(f->coeffs[i], g->coeffs[i], prime);
  }
}

void poly_add(poly *r, poly *f, poly *g) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = f->coeffs[i] + g->coeffs[i];
  }
}

void poly_sub32(poly *r, poly *f, poly *g, int32_t prime) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = sub32(f->coeffs[i], g->coeffs[i], prime);
  }
}

void poly_sub(poly *r, poly *f, poly *g) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = f->coeffs[i] - g->coeffs[i];
  }
}

/// @brief Reduce the input modulo 2**d into signed interval [-2**d, 2**d-1]
void poly_reduce_pow2(poly *r, poly *f, int32_t d) {
  int32_t pow2div2 = (int32_t)1 << (d-1); // 2^d/2
  int32_t pow2 = (int32_t)1 << d;         // 2^d
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = ((f->coeffs[i] + pow2div2) & (pow2-1)) - pow2div2;
  }
}

/// @brief Reduce the input modulo prime into signed interval [-(prime-1)/2, (prime-1)/2-1]
void poly_conditional_final_reduce32(poly *r, int32_t prime) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = conditional_final_reduce32(r->coeffs[i], prime);
  }
}

void poly_round_xtoy(poly *r, const poly *f, int32_t x, int32_t y) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = f->coeffs[i] + ((int32_t)1 << (x-(y+1))); // Add constant x/(2*y)
    r->coeffs[i] >>= (x-y);                                  // Divide by x/y and floor
    r->coeffs[i] &= ((int32_t)1 << y)-1;                     // Reduce mod y
  }
}

void poly_compress(poly *r, int32_t x) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] >>= x;
  }
}

void poly_decompress(poly *r, int32_t x) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] <<= x;
  }
}
