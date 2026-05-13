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

#include "pke.h"

int32_t rrlwr_pke_zetas[RRLWR_N] = RRLWR_KEM_ZETAS;
#ifdef PRECOMPUTE_TWIST
int32_t precomputed_twist[RRLWR_N] = {403905201, -403348161, -66371421, 66928461, 133626020, -133068980, -300247038, 300804078, -117505540, 118062580, -347058431, 347615471, 902218, -345178, 390549999, -389992959, 39812103, -39255063, -96655799, 97212839, -172660914, 173217954, -289881212, 290438252, 199364025, -198806985, 2418007, -1860967, -530070296, 530627336, 519621497, -519064457, 93727338, -93170298, 240048030, -239490990, 73001961, -72444921, -300824911, 301381951, 76681729, -76124689, -227812543, 228369583, 94796637, -94239597, -121921241, 122478281, -22994915, 23551955, -148310860, 148867900, 225965532, -225408492, -134941793, 135498833, 227843057, -227286017, -35874344, 36431384, 324335842, -323778802, -250118201, 250675241, 475909322, -475352282, 363951894, -363394854, 385799325, -385242285, -218094760, 218651800, 198357559, -197800519, 534968254, -534411214, -514371882, 514928922, -369263088, 369820128, 267069560, -266512520, -37236156, 37793196, -175636305, 176193345, 506593685, -506036645, -469797415, 470354455, 204971739, -204414699, 247968200, -247411160, 499704399, -499147359, 216056503, -215499463, -169442834, 169999874, 132993396, -132436356, 106322988, -105765948, 195081012, -194523972, 45073110, -44516070, 330100924, -329543884, 317036508, -316479468, 518816756, -518259716, -15497573, 16054613, -500741754, 501298794, 215992096, -215435056, -382957781, 383514821, 230555203, -229998163, -320273567, 320830607, -182129223, 182686263};
#endif

void poly_subp(poly *r, const poly *f, const poly *g) {
  for(unsigned int i = 0; i < RRLWR_N; i++) {
    int32_t h2 = (int32_t)1 << (RRLWR_PKE_LOGP-(RRLWR_PKE_LOGT+1));      // p/(2*t)
    h2 -= (int32_t)1 << (RRLWR_PKE_LOGQ-(RRLWR_PKE_LOGP+1));             // p/(2*t) - q/(2*p)
    r->coeffs[i] = (f->coeffs[i] - g->coeffs[i] + h2) & (RRLWR_PKE_P-1); // Reduce mod p
  }
}

void poly_add_msg(poly *r, const unsigned char msg[RRLWR_PKE_MESSAGE_LEN]) {
  int32_t mj;

  for(unsigned int i = 0; i < RRLWR_PKE_ELL; i++) {
    for(unsigned int j = 0; j < RRLWR_N; j++) {
      mj = (msg[i*RRLWR_N/8 + (j >> 3)] >> (j & 0x7)) & 1;                 // Select coefficient j of polynomial i
      r[i].coeffs[j] += (int32_t)1 << (RRLWR_PKE_LOGQ-(RRLWR_PKE_LOGP+1)); // Add q/(2*p)
      r[i].coeffs[j] += -mj & ((int32_t)1 << (RRLWR_PKE_LOGP-1));          // Add p/2 if message bit is 1, otherwise not
      r[i].coeffs[j] &= (RRLWR_PKE_P-1);                                   // Reduce mod p
    }
  }
}

uint8_t ct_cmp(const unsigned char *c1, const unsigned char *c2) {
  uint8_t r = 0;
  for(unsigned int i = 0; i < RRLWR_KEM_CT_LEN; i++) {
    r |= (c1[i] ^ c2[i]); // Only non-zero if c1 != c2
  }

  return -((-(uint32_t)r) >> 31);
}

int pke_keygen(unsigned char pk[RRLWR_PKE_PK_LEN], unsigned char sk[RRLWR_PKE_SK_LEN],
               const unsigned char seedA[RRLWR_PKE_SEED_A_LEN], const unsigned char seedS[RRLWR_SEED_S_LEN]) {

  ring_element a, s, b;

  // Generate a with coefficients in [-q/2+1, q/2]
  ring_uniform(&a, RRLWR_PKE_LOGQ, seedA, RRLWR_PKE_SEED_A_LEN);

  // Generate s with coefficients in [-2, 1] and pack
  ring_uniform(&s, RRLWR_PKE_LOG_ETA+1, seedS, RRLWR_SEED_S_LEN);
  ring_pack(sk, &s, RRLWR_PKE_LOG_ETA+1);

  // Compute A*s
  ring_mul32(b.x, &a, &s, RRLWR_K, RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV, RRLWR_NTTINV_FINALCONST, RRLWR_KEM_RMODPRIME, RRLWR_KEM_2RMODPRIME, rrlwr_pke_zetas);

  // Compute b = round(p/q*A*s) computed as (A*s + q/(2*p)) >> (eq - ep)
  ring_round_xtoy(&b, &b, RRLWR_PKE_LOGQ, RRLWR_PKE_LOGP);

  // Pack the public key
  for (unsigned int i = 0; i < RRLWR_PKE_SEED_A_LEN; i++) {
    pk[i] = seedA[i];
  }
  ring_pack(pk + RRLWR_PKE_SEED_A_LEN, &b, RRLWR_PKE_LOGP);

  return 0;
}

int pke_encrypt(unsigned char ct[RRLWR_PKE_CT_LEN], const unsigned char pk[RRLWR_PKE_PK_LEN],
                const unsigned char m[RRLWR_PKE_MESSAGE_LEN], const unsigned char seedSp[RRLWR_SEED_S_LEN]) {

  ring_element a, sp, b, bp;
  poly vp[RRLWR_PKE_ELL];
  const unsigned char *seedA = &pk[0];

  // Generate a with coefficients in [-q/2+1, q/2]
  ring_uniform(&a, RRLWR_PKE_LOGQ, seedA, RRLWR_PKE_SEED_A_LEN);

  // Generate s_prime with coefficients in [-2, 1]
  ring_uniform(&sp, RRLWR_PKE_LOG_ETA+1, seedSp, RRLWR_SEED_S_LEN);

  // Compute A*s_prime
  ring_mul32(bp.x, &a, &sp, RRLWR_K, RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV, RRLWR_NTTINV_FINALCONST, RRLWR_KEM_RMODPRIME, RRLWR_KEM_2RMODPRIME, rrlwr_pke_zetas);

  // Compute b_prime = round(p/q*A*s_prime)
  ring_round_xtoy(&bp, &bp, RRLWR_PKE_LOGQ, RRLWR_PKE_LOGP);

  // Compute v_prime = b*s_prime
  ring_unpack(&b, pk + RRLWR_PKE_SEED_A_LEN, RRLWR_PKE_LOGP);
  ring_ntt32(&b, RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV, rrlwr_pke_zetas);
  ring_mul_invntt32(vp, &b, &sp, RRLWR_PKE_ELL, RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV, RRLWR_NTTINV_FINALCONST, RRLWR_KEM_RMODPRIME, RRLWR_KEM_2RMODPRIME, rrlwr_pke_zetas); // Do not repeat NTT(s)

  // Convert message to polynomial representation and compute (v_prime + q/(2*p) + p/2*m) mod p
  poly_add_msg(vp, m);

  // Round from R_p to R_t with result in [-t/2+1, t/2] and pack into ciphertext buffer
  for(unsigned int i = 0; i < RRLWR_PKE_ELL; i++) {
    poly_compress(&vp[i], RRLWR_PKE_LOGP - RRLWR_PKE_LOGT); // Multiply by t/p and floor
    poly_pack(ct + i*RRLWR_PKE_PACKED_POLYT_LEN, &vp[i], RRLWR_PKE_LOGT);
  }

  // Return ct = (cm = vp, b) and sk = s
  ring_pack(ct + RRLWR_PKE_ELL*RRLWR_PKE_PACKED_POLYT_LEN, &bp, RRLWR_PKE_LOGP);

  return 0;
}

int pke_decrypt(unsigned char m[RRLWR_PKE_MESSAGE_LEN], 
                const unsigned char ct[RRLWR_PKE_CT_LEN], const unsigned char sk[RRLWR_PKE_SK_LEN]) {
  ring_element bp, s; 
  poly v[RRLWR_PKE_ELL];
  poly cm[RRLWR_PKE_ELL];

  // Unpack s, b and cm and decompress cm
  ring_unpack(&s, sk, RRLWR_PKE_LOG_ETA+1);
  ring_unpack(&bp, ct + RRLWR_PKE_ELL*RRLWR_PKE_PACKED_POLYT_LEN, RRLWR_PKE_LOGP);
  for(unsigned int i = 0; i < RRLWR_PKE_ELL; i++) {
    poly_unpack(&cm[i], ct + i*RRLWR_PKE_PACKED_POLYT_LEN, RRLWR_PKE_LOGT);
    poly_decompress(&cm[i], RRLWR_PKE_LOGP - RRLWR_PKE_LOGT); // Multiply by p/t
  }

  // Compute v = b_prime*s
  ring_mul32(v, &bp, &s, RRLWR_PKE_ELL, RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV, RRLWR_NTTINV_FINALCONST, RRLWR_KEM_RMODPRIME, RRLWR_KEM_2RMODPRIME, rrlwr_pke_zetas); 

  for(unsigned int i = 0; i < RRLWR_PKE_ELL; i++) {
    poly_subp(&v[i], &v[i], &cm[i]);                     // Compute (v - (p/t)*cm) mod p
    poly_round_xtoy(&v[i], &v[i], RRLWR_PKE_LOGP, 1);    // Round to a single-bit message polynomial m'
    poly_pack(m + i*RRLWR_PKE_PACKED_POLY1_LEN, &v[i], 1); // Convert the message polynomial to a bit string
  }

  return 0;
}