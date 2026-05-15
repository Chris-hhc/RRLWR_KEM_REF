/*
 * Sliding-window matrix-vector multiplication benchmark.
 *
 * Timed region starts at public A expansion into the sliding-window layout:
 *
 *   NTT path:
 *     seedA -> Awin_NTT, then Awin_NTT * s, including s NTT, matrix-vector
 *     base multiplication, inverse NTT, and final reduction.
 *
 *   coefficient q13 paths:
 *     seedA -> Awin_q13, then Awin_q13 * s using either generic Toom-Cook
 *     polynomial multiplication or the small-secret negative wrapped
 *     convolution path.
 *
 * DRNG and secret sampling are outside the timed region. Rotating A seeds are
 * used so repeated runs do not only exercise one fixed public A.
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
static unsigned char seedA_list[NUMBER_OF_TESTS][RRLWR_PKE_SEED_A_LEN];
static ring_element secret_list[NUMBER_OF_TESTS];

DRNG_ctx drng_algorithm;
extern int32_t rrlwr_pke_zetas[RRLWR_N];

static void init_inputs(void) {
  const unsigned char seed0[32] = {0};
  unsigned char seedS[RRLWR_SEED_S_LEN];

  init_random_number(&drng_algorithm, seed0, sizeof(seed0));

  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    get_random_number(&drng_algorithm, seedA_list[i],
                      8 * RRLWR_PKE_SEED_A_LEN);
    get_random_number(&drng_algorithm, seedS, 8 * RRLWR_SEED_S_LEN);
    ring_uniform(&secret_list[i], RRLWR_PKE_LOG_ETA + 1,
                 seedS, RRLWR_SEED_S_LEN);
  }
}

static void bench_ntt_Awin(const char *name, int ncoeffs) {
  ring_element_Awin_ntt aw;
  ring_element s;
  poly r[RRLWR_K];

  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    s = secret_list[i];
    ring_uniform_Awin_ntt(&aw, RRLWR_PKE_LOGQ,
                          seedA_list[i], RRLWR_PKE_SEED_A_LEN,
                          RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV,
                          RRLWR_KEM_RMODPRIME, RRLWR_KEM_2RMODPRIME,
                          rrlwr_pke_zetas);
    ring_mul32_Awin(r, &aw, &s, ncoeffs,
                    RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV,
                    RRLWR_NTTINV_FINALCONST, rrlwr_pke_zetas);
    t[i] = cpucycles();
  }

  print_results(name, t, NUMBER_OF_TESTS);
}

static void bench_q13_toom_Awin(const char *name, int ncoeffs) {
  ring_element_Awin_q13 aw;
  poly r[RRLWR_K];

  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    ring_uniform_Awin_q13(&aw, RRLWR_PKE_LOGQ,
                          seedA_list[i], RRLWR_PKE_SEED_A_LEN);
    ring_mul_q13_Awin(r, &aw, &secret_list[i], ncoeffs);
    t[i] = cpucycles();
  }

  print_results(name, t, NUMBER_OF_TESTS);
}

static void bench_q13_smallsecret_Awin(const char *name, int ncoeffs) {
  ring_element_Awin_q13 aw;
  poly r[RRLWR_K];

  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    ring_uniform_Awin_q13(&aw, RRLWR_PKE_LOGQ,
                          seedA_list[i], RRLWR_PKE_SEED_A_LEN);
    ring_mul_q13_Awin_smallsecret(r, &aw, &secret_list[i], ncoeffs);
    t[i] = cpucycles();
  }

  print_results(name, t, NUMBER_OF_TESTS);
}

static void bench_q13_smallsecret_ct_Awin(const char *name, int ncoeffs) {
  ring_element_Awin_q13 aw;
  poly r[RRLWR_K];

  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    ring_uniform_Awin_q13(&aw, RRLWR_PKE_LOGQ,
                          seedA_list[i], RRLWR_PKE_SEED_A_LEN);
    ring_mul_q13_Awin_smallsecret_ct(r, &aw, &secret_list[i], ncoeffs);
    t[i] = cpucycles();
  }

  print_results(name, t, NUMBER_OF_TESTS);
}

int main(void) {
  init_inputs();

  bench_ntt_Awin("Awin NTT sample + mv + intt/finalreduce, full K", RRLWR_K);
  bench_q13_toom_Awin("Awin q13 sample + Toom-Cook mv, full K", RRLWR_K);
  bench_q13_smallsecret_Awin("Awin q13 sample + smallsecret negwrap mv, full K", RRLWR_K);
  bench_q13_smallsecret_ct_Awin("Awin q13 sample + CT smallsecret negwrap mv, full K", RRLWR_K);

  bench_ntt_Awin("Awin NTT sample + mv + intt/finalreduce, ell", RRLWR_PKE_ELL);
  bench_q13_toom_Awin("Awin q13 sample + Toom-Cook mv, ell", RRLWR_PKE_ELL);
  bench_q13_smallsecret_Awin("Awin q13 sample + smallsecret negwrap mv, ell", RRLWR_PKE_ELL);
  bench_q13_smallsecret_ct_Awin("Awin q13 sample + CT smallsecret negwrap mv, ell", RRLWR_PKE_ELL);

  return 0;
}
