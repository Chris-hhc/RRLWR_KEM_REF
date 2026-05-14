/*
 * Benchmark q13 Awin small-secret multiplication paths.
 *
 * This compares the existing sparse-index non-CT path against the CT path for:
 *   - full RRLWR_K output rows
 *   - RRLWR_PKE_ELL output rows
 *   - fixed secret
 *   - rotating pre-generated secrets
 */

#include <stdint.h>

#include "cpucycles.h"
#include "drng.h"
#include "parameters.h"
#include "ring.h"
#include "speed_print.h"
#include "uniform.h"

#define NUMBER_OF_TESTS 1000

static uint64_t t[NUMBER_OF_TESTS];
static ring_element secrets[NUMBER_OF_TESTS];

DRNG_ctx drng_algorithm;

static void bench_fixed_secret(const char *name,
                               int use_ct,
                               poly *r,
                               const ring_element_Awin_q13 *aw,
                               const ring_element *s,
                               int ncoeffs) {
  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    if(use_ct) {
      ring_mul_q13_Awin_smallsecret_ct(r, aw, s, ncoeffs);
    } else {
      ring_mul_q13_Awin_smallsecret(r, aw, s, ncoeffs);
    }
    t[i] = cpucycles();
  }

  print_results(name, t, NUMBER_OF_TESTS);
}

static void bench_rotating_secret(const char *name,
                                  int use_ct,
                                  poly *r,
                                  const ring_element_Awin_q13 *aw,
                                  int ncoeffs) {
  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    if(use_ct) {
      ring_mul_q13_Awin_smallsecret_ct(r, aw, &secrets[i], ncoeffs);
    } else {
      ring_mul_q13_Awin_smallsecret(r, aw, &secrets[i], ncoeffs);
    }
    t[i] = cpucycles();
  }

  print_results(name, t, NUMBER_OF_TESTS);
}

int main(void) {
  const unsigned char seed0[32] = {0};
  unsigned char seedA[RRLWR_PKE_SEED_A_LEN];
  unsigned char seedS[RRLWR_SEED_S_LEN];
  ring_element_Awin_q13 aw;
  poly r[RRLWR_K];

  init_random_number(&drng_algorithm, seed0, sizeof(seed0));
  get_random_number(&drng_algorithm, seedA, 8 * RRLWR_PKE_SEED_A_LEN);
  ring_uniform_Awin_q13(&aw, RRLWR_PKE_LOGQ, seedA, RRLWR_PKE_SEED_A_LEN);

  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    get_random_number(&drng_algorithm, seedS, 8 * RRLWR_SEED_S_LEN);
    ring_uniform(&secrets[i], RRLWR_PKE_LOG_ETA + 1, seedS, RRLWR_SEED_S_LEN);
  }

  bench_fixed_secret("non-CT smallsecret full K, fixed s",
                     0, r, &aw, &secrets[0], RRLWR_K);
  bench_fixed_secret("CT smallsecret full K, fixed s",
                     1, r, &aw, &secrets[0], RRLWR_K);

  bench_fixed_secret("non-CT smallsecret ell, fixed s",
                     0, r, &aw, &secrets[0], RRLWR_PKE_ELL);
  bench_fixed_secret("CT smallsecret ell, fixed s",
                     1, r, &aw, &secrets[0], RRLWR_PKE_ELL);

  bench_rotating_secret("non-CT smallsecret full K, rotating s",
                        0, r, &aw, RRLWR_K);
  bench_rotating_secret("CT smallsecret full K, rotating s",
                        1, r, &aw, RRLWR_K);

  bench_rotating_secret("non-CT smallsecret ell, rotating s",
                        0, r, &aw, RRLWR_PKE_ELL);
  bench_rotating_secret("CT smallsecret ell, rotating s",
                        1, r, &aw, RRLWR_PKE_ELL);

  return 0;
}
