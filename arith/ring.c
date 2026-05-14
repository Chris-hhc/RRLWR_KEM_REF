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

#include "ring.h"
#include "uniform.h"

#ifndef NDEBUG
#include <stdio.h>
#include <stdlib.h>
#endif

/// @brief Compute NTT(y+2) in Montgomery domain, on-the-fly only if PRECOMPUTE_TWIST is not defined
static void compute_yp2(poly *r, int32_t prime, int32_t primeinv, int32_t fp_zetas[RRLWR_N], int32_t oneR, int32_t twoR)
{
  #ifdef PRECOMPUTE_TWIST
    (void)prime;
    (void)primeinv;
    (void)fp_zetas;
    (void)oneR;
    (void)twoR;
    #ifdef CRT
    if(prime == RRLWR_SIGN_PRIME1) {
      for(unsigned int i = 0; i < RRLWR_N; i++) {
        r->coeffs[i] = precomputed_twist1[i];
      }
    } else if (prime == RRLWR_SIGN_PRIME2) {
      for(unsigned int i = 0; i < RRLWR_N; i++) {
        r->coeffs[i] = precomputed_twist2[i];
      }
    }
    #else
    for(unsigned int i = 0; i < RRLWR_N; i++) {
      r->coeffs[i] = precomputed_twist[i];
    }
    #endif
  #else
    // Initialize the polynomial y+2 in Montgomery domain
    r->coeffs[0] = twoR; // 2
    r->coeffs[1] = oneR; // 1
    for(unsigned int i = 2; i < RRLWR_N; i++) {
      r->coeffs[i] = 0;
    }
    poly_ntt32(r, prime, primeinv, fp_zetas);
  #endif
}

void ring_ntt32(ring_element *r, int32_t prime, int32_t primeinv, int32_t fp_zetas[RRLWR_N]) {
  for(int i = 0; i < RRLWR_K; i++) {
    poly_ntt32(&r->x[i], prime, primeinv, fp_zetas);
  }
}

void ring_uniform_Awin_ntt(ring_element_Awin_ntt *aw,
                           int32_t bitlen,
                           const unsigned char *seed,
                           int32_t seed_len,
                           int32_t prime,
                           int32_t primeinv,
                           int32_t oneR,
                           int32_t twoR,
                           int32_t fp_zetas[RRLWR_N]) {
  poly a, yp2;

  if(RRLWR_K > 1) {
    compute_yp2(&yp2, prime, primeinv, fp_zetas, oneR, twoR);
  }

  for(unsigned char u = 0; u < RRLWR_K; u++) {
    poly_uniform(&a, bitlen, seed, seed_len, u);
    poly_ntt32(&a, prime, primeinv, fp_zetas);

    aw->x[RRLWR_K - 1 - u] = a;

    if(u > 0) {
      poly_basemul32(&aw->x[2 * RRLWR_K - 1 - u], &a, &yp2, prime, primeinv);
    }
  }
}

void ring_to_Awin_q13(ring_element_Awin_q13 *aw, const ring_element *a) {
  for(unsigned char u = 0; u < RRLWR_K; u++) {
    aw->x[RRLWR_K - 1 - u] = a->x[u];
  }

  for(unsigned char u = 1; u < RRLWR_K; u++) {
    int base = RRLWR_K - 1 - u;
    int dst = 2 * RRLWR_K - 1 - u;

    poly_mul_yplus2_q13(&aw->x[dst], &aw->x[base]);
  }
}

void ring_uniform_Awin_q13(ring_element_Awin_q13 *aw,
                           int32_t bitlen,
                           const unsigned char *seed,
                           int32_t seed_len) {
  poly a;

  for(unsigned char u = 0; u < RRLWR_K; u++) {
    poly_uniform(&a, bitlen, seed, seed_len, u);
    aw->x[RRLWR_K - 1 - u] = a;
  }

  for(unsigned char u = 1; u < RRLWR_K; u++) {
    int base = RRLWR_K - 1 - u;
    int dst = 2 * RRLWR_K - 1 - u;

    poly_mul_yplus2_q13(&aw->x[dst], &aw->x[base]);
  }
}

/// @brief Full ring multiplication assuming its inputs are already in NTT domain.
///        The output element r is not in NTT domain and is fully reduced with all coefficients in [-q/2, q/2+1]
// void ring_mul_invntt32(poly *r, ring_element *a, ring_element *b, int ncoeffs, int32_t prime, int32_t primeinv, int32_t finalconst, int32_t oneR, int32_t twoR, int32_t fp_zetas[RRLWR_N]) {
//   poly yp2, t;

//   if(RRLWR_K > 1) { // Not required if K = 1
//     compute_yp2(&yp2, prime, primeinv, fp_zetas, oneR, twoR); // NTT(y+2) in Montgomery domain
//   }

//   int rindex = ncoeffs-1;

//   // Remaining rows of matrix multiplication that are needed according to ncoeffs
//   for(int i = RRLWR_K-1; i > RRLWR_K-1-ncoeffs; i--) {

//     for(int j = 0; j < RRLWR_N; j++) {
//       (r+rindex)->coeffs[j] = 0;
//     }

//     // Process row elements not multiplied by (y+2)
//     for(int j = 0; j < i+1; j++) {
//       poly_basemul32(&t, &a->x[i-j], &b->x[j], prime, primeinv); // Introduces a Montgomery factor R^-1
//       poly_add32(&r[rindex], &r[rindex], &t, prime);
//     }

//     // Multiply matrix element by (y+2)
//     if (i+1 < RRLWR_K) {
//       poly_basemul32(&a->x[i+1], &a->x[i+1], &yp2, prime, primeinv); // Montgomery factor R cancelled from yp2
//     }

//     // Process row elements multiplied by (y+2)
//     for(int j = i+1; j < RRLWR_K; j++) {
//       poly_basemul32(&t, &a->x[i-j+RRLWR_K], &b->x[j], prime, primeinv); // Introduces a Montgomery factor R^-1
//       poly_add32(&r[rindex], &r[rindex], &t, prime);
//     }

//     rindex--;
//   }

//   // Transfrom to polynomial domain
//   for(int i = 0; i < ncoeffs; i++) {
//     poly_invntt32(&r[i], prime, primeinv, finalconst, fp_zetas); // Removes the factor R^-1 in the final multiplication
//     poly_conditional_final_reduce32(&r[i], prime); // Reduce to unique representation in [-(p-1)/2+1, (p-1)/2]
//   }
// }
void ring_mul_invntt32(poly *r,
                       ring_element *a,
                       ring_element *b,
                       int ncoeffs,
                       int32_t prime,
                       int32_t primeinv,
                       int32_t finalconst,
                       int32_t oneR,
                       int32_t twoR,
                       int32_t fp_zetas[RRLWR_N]) {
  poly yp2, t;

  /*
   * Aline =
   *   [a[k-1], a[k-2], ..., a[0],
   *    theta*a[k-1], theta*a[k-2], ..., theta*a[1]]
   *
   * Only the top ncoeffs rows are needed, so we only need theta*a[k-1]
   * down to theta*a[k-ncoeffs+1].
   *
   * The array is allocated at max size for C89/C90 friendliness.
   * Actually initialized range is [0, RRLWR_K + ncoeffs - 2].
   */
  poly Aline[2 * RRLWR_K - 1];

  /*
   * First half: reversed a.
   *
   * Aline[0]       = a[k-1]
   * Aline[1]       = a[k-2]
   * ...
   * Aline[k-1]     = a[0]
   */
  for (int u = 0; u < RRLWR_K; u++) {
    Aline[u] = a->x[RRLWR_K - 1 - u];
  }

  /*
   * Second half: theta * a, also in reversed order.
   *
   * Aline[k]       = theta * a[k-1]
   * Aline[k+1]     = theta * a[k-2]
   * ...
   *
   * For top ncoeffs rows, only ncoeffs-1 wrapped entries can appear.
   */
  if (RRLWR_K > 1 && ncoeffs > 1) {
    compute_yp2(&yp2, prime, primeinv, fp_zetas, oneR, twoR);

    for (int u = 1; u < ncoeffs; u++) {
      int src = RRLWR_K - u;
      int dst = RRLWR_K - 1 + u;

      /*
       * Same Montgomery-factor convention as the original code:
       * original did:
       *   a->x[src] = a->x[src] * yp2
       * before using it in wrapped row entries.
       */
      poly_basemul32(&Aline[dst], &a->x[src], &yp2, prime, primeinv);
    }
  }

  int row_min = RRLWR_K - ncoeffs;

  /*
   * Compute rows i = k-1, ..., k-ncoeffs.
   *
   * Preserve the original output order:
   *   row k-1         -> r[ncoeffs-1]
   *   row k-ncoeffs   -> r[0]
   */
  for (int i = RRLWR_K - 1; i >= row_min; i--) {
    int out = i - row_min;

    /*
     * row points to:
     *   [a[i], a[i-1], ..., a[0],
     *    theta*a[k-1], theta*a[k-2], ...]
     */
    poly *row = &Aline[RRLWR_K - 1 - i];

    poly_basemul32(&r[out], &row[0], &b->x[0], prime, primeinv);

    for (int j = 1; j < RRLWR_K; j++) {
      poly_basemul32(&t, &row[j], &b->x[j], prime, primeinv);
      poly_add32(&r[out], &r[out], &t, prime);
    }
  }

  /*
   * Transform requested output coefficients back to polynomial domain.
   */
  for (int i = 0; i < ncoeffs; i++) {
    poly_invntt32(&r[i], prime, primeinv, finalconst, fp_zetas);
    poly_conditional_final_reduce32(&r[i], prime);
  }
}

void ring_mul_invntt32_Awin(poly *r,
                            const ring_element_Awin_ntt *a,
                            ring_element *b,
                            int ncoeffs,
                            int32_t prime,
                            int32_t primeinv,
                            int32_t finalconst,
                            int32_t fp_zetas[RRLWR_N]) {
  poly t;
  int row_min = RRLWR_K - ncoeffs;

  for (int i = RRLWR_K - 1; i >= row_min; i--) {
    int out = i - row_min;
    const poly *row = &a->x[RRLWR_K - 1 - i];

    poly_basemul32(&r[out], &row[0], &b->x[0], prime, primeinv);

    for (int j = 1; j < RRLWR_K; j++) {
      poly_basemul32(&t, &row[j], &b->x[j], prime, primeinv);
      poly_add32(&r[out], &r[out], &t, prime);
    }
  }

  for (int i = 0; i < ncoeffs; i++) {
    poly_invntt32(&r[i], prime, primeinv, finalconst, fp_zetas);
    poly_conditional_final_reduce32(&r[i], prime);
  }
}

// void ring_mul_invntt32(poly *r,
//                        ring_element *a,
//                        ring_element *b,
//                        int ncoeffs,
//                        int32_t prime,
//                        int32_t primeinv,
//                        int32_t finalconst,
//                        int32_t oneR,
//                        int32_t twoR,
//                        int32_t fp_zetas[RRLWR_N]) {
//   poly yp2, t;
//   poly a_yp2[RRLWR_K];

//   // if (RRLWR_K > 1) { // Not required if K = 1
//   //   compute_yp2(&yp2, prime, primeinv, fp_zetas, oneR, twoR); // NTT(y+2) in Montgomery domain

//   //   // Precompute a * (y+2) outside the main loop.
//   //   // a_yp2[0] is not needed for wrapped terms, but computing from 1 is enough.
//   //   for (int i = 1; i < RRLWR_K; i++) {
//   //     poly_basemul32(&a_yp2[i], &a->x[i], &yp2, prime, primeinv);
//   //     // Montgomery factor R cancelled from yp2, same as the original in-loop update.
//   //   }
//   // }

//   int rindex = ncoeffs - 1;

//   // Remaining rows of matrix multiplication that are needed according to ncoeffs
//   for (int i = RRLWR_K - 1; i > RRLWR_K - 1 - ncoeffs; i--) {

//     // for (int j = 0; j < RRLWR_N; j++) {
//     //   r[rindex].coeffs[j] = 0;
//     // }

//     // Process row elements not multiplied by (y+2)
//     for (int j = 0; j < i + 1; j++) {
//       poly_basemul32(&t, &a->x[i - j], &b->x[j], prime, primeinv);
//       // Introduces a Montgomery factor R^-1
//       poly_add32(&r[rindex], &r[rindex], &t, prime);
//     }

//     // Process row elements multiplied by (y+2)
//     // Use precomputed a_yp2 instead of modifying a->x in place.
//     for (int j = i + 1; j < RRLWR_K; j++) {
//       poly_basemul32(&t, &a_yp2[i - j + RRLWR_K], &b->x[j], prime, primeinv);
//       // Introduces a Montgomery factor R^-1
//       poly_add32(&r[rindex], &r[rindex], &t, prime);
//     }

//     rindex--;
//   }

//   // Transform to polynomial domain
//   for (int i = 0; i < ncoeffs; i++) {
//     // poly_invntt32(&r[i], prime, primeinv, finalconst, fp_zetas);
//     // Removes the factor R^-1 in the final multiplication
//     // poly_conditional_final_reduce32(&r[i], prime);
//     // Reduce to unique representation in [-(p-1)/2+1, (p-1)/2]
//   }
// }

/// @brief Full ring multiplication assuming its inputs are not yet in NTT domain.
///        The output element r is not in NTT domain and is fully reduced with all coefficients in [-q/2, q/2+1]
void ring_mul32(poly *r, ring_element *a, ring_element *b, int ncoeffs, int32_t prime, int32_t primeinv, int32_t finalconst, int32_t oneR, int32_t twoR, int32_t fp_zetas[RRLWR_N])
{
  ring_ntt32(a, prime, primeinv, fp_zetas);
  ring_ntt32(b, prime, primeinv, fp_zetas);
  ring_mul_invntt32(r, a, b, ncoeffs, prime, primeinv, finalconst, oneR, twoR, fp_zetas);
}

void ring_mul32_Awin(poly *r, const ring_element_Awin_ntt *a, ring_element *b, int ncoeffs, int32_t prime, int32_t primeinv, int32_t finalconst, int32_t fp_zetas[RRLWR_N])
{
  ring_ntt32(b, prime, primeinv, fp_zetas);
  ring_mul_invntt32_Awin(r, a, b, ncoeffs, prime, primeinv, finalconst, fp_zetas);
}

void ring_mul_q13_Awin(poly *r,
                       const ring_element_Awin_q13 *a,
                       const ring_element *b,
                       int ncoeffs) {
  int row_min = RRLWR_K - ncoeffs;

  for(int i = RRLWR_K - 1; i >= row_min; i--) {
    int out = i - row_min;
    const poly *row = &a->x[RRLWR_K - 1 - i];

    poly_mul_q13(&r[out], &row[0], &b->x[0]);

    for(int j = 1; j < RRLWR_K; j++) {
      poly_macc_q13(&r[out], &row[j], &b->x[j]);
    }
  }
}

void ring_mul_q13(poly *r,
                  const ring_element *a,
                  const ring_element *b,
                  int ncoeffs) {
  ring_element_Awin_q13 aw;

  ring_to_Awin_q13(&aw, a);
  ring_mul_q13_Awin(r, &aw, b, ncoeffs);
}

static void ring_row_macc_q13_Awin_smallsecret(poly *out,
                                               const int32_t (*row)[RRLWR_N],
                                               const uint16_t idx_p1[RRLWR_K][RRLWR_N],
                                               const uint16_t idx_m1[RRLWR_K][RRLWR_N],
                                               const uint16_t idx_p2[RRLWR_K][RRLWR_N],
                                               const uint16_t idx_m2[RRLWR_K][RRLWR_N],
                                               const unsigned int count_p1[RRLWR_K],
                                               const unsigned int count_m1[RRLWR_K],
                                               const unsigned int count_p2[RRLWR_K],
                                               const unsigned int count_m2[RRLWR_K]) {
  int32_t acc[RRLWR_N];

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    acc[i] = 0;
  }

  /*
   * This ref fast path branches on small secret coefficients; replace with
   * constant-time masks if required.
   */
  for(unsigned int j = 0; j < RRLWR_K; j++) {
    for(unsigned int l = 0; l < count_p1[j]; l++) {
      unsigned int c = idx_p1[j][l];

      for(unsigned int t = 0; t < RRLWR_N - c; t++) {
        acc[t + c] += row[j][t];
      }
      for(unsigned int t = RRLWR_N - c; t < RRLWR_N; t++) {
        acc[t + c - RRLWR_N] -= row[j][t];
      }
    }

    for(unsigned int l = 0; l < count_m1[j]; l++) {
      unsigned int c = idx_m1[j][l];

      for(unsigned int t = 0; t < RRLWR_N - c; t++) {
        acc[t + c] -= row[j][t];
      }
      for(unsigned int t = RRLWR_N - c; t < RRLWR_N; t++) {
        acc[t + c - RRLWR_N] += row[j][t];
      }
    }

    for(unsigned int l = 0; l < count_p2[j]; l++) {
      unsigned int c = idx_p2[j][l];

      for(unsigned int t = 0; t < RRLWR_N - c; t++) {
        acc[t + c] += row[j][t] << 1;
      }
      for(unsigned int t = RRLWR_N - c; t < RRLWR_N; t++) {
        acc[t + c - RRLWR_N] -= row[j][t] << 1;
      }
    }

    for(unsigned int l = 0; l < count_m2[j]; l++) {
      unsigned int c = idx_m2[j][l];

      for(unsigned int t = 0; t < RRLWR_N - c; t++) {
        acc[t + c] -= row[j][t] << 1;
      }
      for(unsigned int t = RRLWR_N - c; t < RRLWR_N; t++) {
        acc[t + c - RRLWR_N] += row[j][t] << 1;
      }
    }
  }

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    out->coeffs[i] = q13_reduce_u_i64(acc[i]);
  }
}

void ring_mul_q13_Awin_smallsecret(poly *r,
                                   const ring_element_Awin_q13 *a,
                                   const ring_element *s,
                                   int ncoeffs) {
  int32_t awc[RRLWR_K + ncoeffs - 1][RRLWR_N];
  uint16_t idx_p1[RRLWR_K][RRLWR_N];
  uint16_t idx_m1[RRLWR_K][RRLWR_N];
  uint16_t idx_p2[RRLWR_K][RRLWR_N];
  uint16_t idx_m2[RRLWR_K][RRLWR_N];
  unsigned int count_p1[RRLWR_K];
  unsigned int count_m1[RRLWR_K];
  unsigned int count_p2[RRLWR_K];
  unsigned int count_m2[RRLWR_K];
  int row_min = RRLWR_K - ncoeffs;

  for(unsigned int w = 0; w < RRLWR_K + (unsigned int)ncoeffs - 1; w++) {
    for(unsigned int c = 0; c < RRLWR_N; c++) {
      awc[w][c] = q13_center_i32(a->x[w].coeffs[c]);
    }
  }

  for(unsigned int j = 0; j < RRLWR_K; j++) {
    count_p1[j] = 0;
    count_m1[j] = 0;
    count_p2[j] = 0;
    count_m2[j] = 0;

    for(unsigned int c = 0; c < RRLWR_N; c++) {
      int32_t sc = q13_small_secret_i32(s->x[j].coeffs[c]);

#ifndef NDEBUG
      if(!(sc == -2 || sc == -1 || sc == 0 || sc == 1 || sc == 2)) {
        fprintf(stderr, "ring_mul_q13_Awin_smallsecret: non-small coeff %d\n", sc);
        abort();
      }
#endif

      if(sc == 1) {
        idx_p1[j][count_p1[j]++] = (uint16_t)c;
      } else if(sc == -1) {
        idx_m1[j][count_m1[j]++] = (uint16_t)c;
      } else if(sc == 2) {
        idx_p2[j][count_p2[j]++] = (uint16_t)c;
      } else if(sc == -2) {
        idx_m2[j][count_m2[j]++] = (uint16_t)c;
      }
    }
  }

  for(int i = RRLWR_K - 1; i >= row_min; i--) {
    int out = i - row_min;
    const int32_t (*row)[RRLWR_N] = &awc[RRLWR_K - 1 - i];

    ring_row_macc_q13_Awin_smallsecret(&r[out], row,
                                       idx_p1, idx_m1, idx_p2, idx_m2,
                                       count_p1, count_m1, count_p2, count_m2);
  }
}

static inline int32_t ct_mask_eq_i32(int32_t x, int32_t y) {
  uint32_t v = (uint32_t)(x ^ y);

  /*
   * Return 0xffffffff if x == y, otherwise 0x00000000.
   */
  v = (v | (uint32_t)(0u - v)) >> 31;
  return -(int32_t)(v ^ 1u);
}

/*
 * CT shifted MAC for sc in {-2,-1,0,1,2}.
 *
 * The masks depend on sc, but they are computed once per secret coefficient c,
 * not once per dense coefficient t.
 */
static inline void acc_shift_smallsecret_ct_scalar(int32_t acc[RRLWR_N],
                                                   const int32_t dense[RRLWR_N],
                                                   unsigned int c,
                                                   int32_t sc) {
  unsigned int split = RRLWR_N - c;

  int32_t m_p1 = ct_mask_eq_i32(sc,  1);
  int32_t m_m1 = ct_mask_eq_i32(sc, -1);
  int32_t m_p2 = ct_mask_eq_i32(sc,  2);
  int32_t m_m2 = ct_mask_eq_i32(sc, -2);

  for(unsigned int t = 0; t < split; t++) {
    int32_t d = dense[t];
    int32_t d2 = d + d;

    int32_t term = 0;
    term += d  & m_p1;
    term -= d  & m_m1;
    term += d2 & m_p2;
    term -= d2 & m_m2;

    acc[t + c] += term;
  }

  for(unsigned int t = split; t < RRLWR_N; t++) {
    int32_t d = dense[t];
    int32_t d2 = d + d;

    int32_t term = 0;
    term += d  & m_p1;
    term -= d  & m_m1;
    term += d2 & m_p2;
    term -= d2 & m_m2;

    acc[t + c - RRLWR_N] -= term;
  }
}
static void ring_row_macc_q13_Awin_smallsecret_ct(poly *out,
                                                  const int32_t (*row)[RRLWR_N],
                                                  const int32_t swc[RRLWR_K][RRLWR_N]) {
  int32_t acc[RRLWR_N];

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    acc[i] = 0;
  }

  /*
   * Fixed loop structure:
   *   j = 0..K-1
   *   c = 0..N-1
   *
   * No skip for sc == 0 and no branch on sc.
   */
  for(unsigned int j = 0; j < RRLWR_K; j++) {
    for(unsigned int c = 0; c < RRLWR_N; c++) {
      acc_shift_smallsecret_ct_scalar(acc, row[j], c, swc[j][c]);
    }
  }

  for(unsigned int i = 0; i < RRLWR_N; i++) {
    out->coeffs[i] = q13_reduce_u_i32(acc[i]);
  }
}

void ring_mul_q13_Awin_smallsecret_ct(poly *r,
                                      const ring_element_Awin_q13 *a,
                                      const ring_element *s,
                                      int ncoeffs) {
  int32_t awc[2 * RRLWR_K - 1][RRLWR_N];
  int32_t swc[RRLWR_K][RRLWR_N];
  int row_min = RRLWR_K - ncoeffs;

  /*
   * Center Awin coefficients once.
   * Awin is public-derived or dense public data; layout remains sliding-window:
   *   [a[k-1], ..., a[0], theta*a[k-1], ..., theta*a[1]]
   */
  for(unsigned int w = 0; w < 2 * RRLWR_K - 1; w++) {
    for(unsigned int c = 0; c < RRLWR_N; c++) {
      awc[w][c] = q13_center_i32_ct(a->x[w].coeffs[c]);
    }
  }

  /*
   * Center secret coefficients once.
   * Expected values after centering: {-2, -1, 0, 1}; +2 is also supported
   * to match the existing non-CT smallsecret interface.
   */
  for(unsigned int j = 0; j < RRLWR_K; j++) {
    for(unsigned int c = 0; c < RRLWR_N; c++) {
      swc[j][c] = q13_center_i32_ct(s->x[j].coeffs[c]);
    }
  }

  /*
   * Preserve original output order:
   *   row k-1       -> r[ncoeffs-1]
   *   row k-ncoeffs -> r[0]
   */
  for(int i = RRLWR_K - 1; i >= row_min; i--) {
    int out = i - row_min;
    const int32_t (*row)[RRLWR_N] = &awc[RRLWR_K - 1 - i];

    ring_row_macc_q13_Awin_smallsecret_ct(&r[out], row, swc);
  }
}

void ring_round_xtoy(ring_element *r, const ring_element *f, int32_t x, int32_t y) {
  for(unsigned int i = 0; i < RRLWR_K; i++) {
    poly_round_xtoy(&r->x[i], &f->x[i], x, y);
  }
}
