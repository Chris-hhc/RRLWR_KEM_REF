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

#include "packing.h"

/// @brief Pack a polynomial with coefficients in [-bitlen/2, bitlen/2-1] into a byte string of bitlen bits per coefficient
void poly_pack(unsigned char *b, poly *r, int32_t bitlen) {
  unsigned int i;
  unsigned int acc_shift = 0; 
  unsigned int bpos = 0;
  poly rp;
  int32_t *rc = rp.coeffs;
  uint32_t acc = 0;

  // Make all coefficients from r positive
  for(i = 0; i < RRLWR_N; i++) {
    rp.coeffs[i] = (((int32_t)1<<(bitlen-1))-1)-r->coeffs[i]; // Subtract from 2^bitlen/2-1 to move to interval [0, bitlen-1]
    rp.coeffs[i] &= ((int32_t)1<<bitlen)-1; // Remove any remaining sign bits
  }

  while (rc < rp.coeffs + RRLWR_N) {

      // Main packing loop
      for(i = 0; i < 32; i++) {
        acc |= (uint32_t)(rc[i] << acc_shift); // Take the next bitlen bits
        acc_shift += bitlen;
        while (acc_shift >= 8) {
          b[bpos++] = (unsigned char)(acc & 0xFF);
          acc >>= 8;
          acc_shift -= 8;
        }
      }

    rc += 32; // Unpack 32 coefficients at a time
  }
}

void ring_pack(unsigned char *b, ring_element *r, int32_t bitlen) {
  unsigned int offset = bitlen*(RRLWR_N>>3);
  for (unsigned int i = 0; i < RRLWR_K; i++) {
    poly_pack(b+i*offset, &r->x[i], bitlen);
  }
}

/// @brief Unpack a polynomial with coefficients in [-bitlen/2, bitlen/2-1] from a byte string of bitlen bits per coefficient
void poly_unpack(poly *r, const unsigned char *b, int32_t bitlen) {
  unsigned int i;
  int32_t acc_shift = 0; 
  unsigned int bpos = 0;
  int32_t *rc = r->coeffs;
  uint32_t acc = 0;

  while (rc < r->coeffs + RRLWR_N) {

      // Main unpacking loop
      for(i = 0; i < 32; i++) {
        while (acc_shift < bitlen) {
          acc |= ((uint32_t)b[bpos++]) << acc_shift;
          acc_shift += 8;
        }
        rc[i] = (int32_t)acc & (((int32_t)1 << bitlen)-1);
        acc >>= bitlen;
        acc_shift -= bitlen;
      }

    rc += 32; // Unpack 32 coefficients at a time
  }

  // Make all coefficients from r signed
  for(i = 0; i < RRLWR_N; i++) {
    r->coeffs[i] = (((int32_t)1<<(bitlen-1))-1)-r->coeffs[i]; // Subtract from 2^bitlen/2-1 to move to interval [-bitlen/2, bitlen/2-1]
  }
}

void ring_unpack(ring_element *r, const unsigned char *b, int32_t bitlen) {
  unsigned int offset = bitlen*(RRLWR_N>>3);
  for (unsigned int i = 0; i < RRLWR_K; i++) {
    poly_unpack(&r->x[i], b+i*offset, bitlen);
  }
}