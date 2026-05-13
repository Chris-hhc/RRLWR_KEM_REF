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

void ring_round_xtoy(ring_element *r, const ring_element *f, int32_t x, int32_t y) {
  for(unsigned int i = 0; i < RRLWR_K; i++) {
    poly_round_xtoy(&r->x[i], &f->x[i], x, y);
  }
}
