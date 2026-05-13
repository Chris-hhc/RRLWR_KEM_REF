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

#include "uniform.h"

/// @brief Generate a polynomial with coefficients pseudo-randomly generated in the uniform distribution [-bitlen/2, bitlen/2-1]
///        The seed (of length seed_len) is appended with coeff that is used as ring coefficient index.
#define RRLWR_MAX_SAMPLING_BITLEN (24)  // Support sampling bit lengths up to 24 bits, only required to define buffer size
#define RRLWR_MAX_SEED_LEN        (128) // Support seed lengths up to 128 bytes, only required to define buffer size
void poly_uniform(poly *r, int32_t bitlen, const unsigned char *seed, int32_t seed_len, unsigned char coeff) {

  unsigned char xof_bytes_buffer[RRLWR_MAX_SAMPLING_BITLEN*(RRLWR_N >> 3)]; // Allocate maximum length
  unsigned char seed_buffer[RRLWR_MAX_SEED_LEN+1];                          // Allocate maximum length

  // Concatenate the seed with the ring element index i
  for(int32_t i = 0; i < seed_len; i++) {
    seed_buffer[i] = seed[i];
  }
  seed_buffer[seed_len] = coeff;

  // Generate the output byte stream
  RRLWR_XOF(xof_bytes_buffer, bitlen*(RRLWR_N >> 3), seed_buffer, seed_len+1);

  // Unpack into polynomial coefficients
  poly_unpack(r, xof_bytes_buffer, bitlen);
}

/// @brief Generate a ring element with coefficients pseudo-randomly generated in the uniform distribution [-bitlen/2, bitlen/2-1]
void ring_uniform(ring_element *r, int32_t bitlen, const unsigned char *seed, int32_t seed_len) {
  for(unsigned char i = 0; i < RRLWR_K; i++) {
    poly_uniform(&r->x[i], bitlen, seed, seed_len, i);
  }
}