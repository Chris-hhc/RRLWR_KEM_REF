/*
 * Polynomial multiplication benchmark.
 *
 * Compares:
 *   1. NTT(a), NTT(b), pointwise Montgomery multiplication, inverse NTT,
 *      final reduction.
 *   2. q13 coefficient-domain Toom-Cook multiplication with negacyclic fold.
 *
 * Inputs are generated from the current parameters.h sampling routines before
 * the timed region.
 */

#include <stdint.h>

#include "cpucycles.h"
#include "drng.h"
#include "parameters.h"
#include "poly.h"
#include "speed_print.h"
#include "uniform.h"

#define NUMBER_OF_TESTS 1000

static uint64_t t[NUMBER_OF_TESTS];
static poly a_list[NUMBER_OF_TESTS];
static poly b_list[NUMBER_OF_TESTS];

DRNG_ctx drng_algorithm;
extern int32_t rrlwr_pke_zetas[RRLWR_N];

static void init_inputs(void) {
  const unsigned char seed0[32] = {0};
  unsigned char seedA[RRLWR_PKE_SEED_A_LEN];
  unsigned char seedB[RRLWR_PKE_SEED_A_LEN];

  init_random_number(&drng_algorithm, seed0, sizeof(seed0));

  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    get_random_number(&drng_algorithm, seedA, 8 * RRLWR_PKE_SEED_A_LEN);
    get_random_number(&drng_algorithm, seedB, 8 * RRLWR_PKE_SEED_A_LEN);

    poly_uniform(&a_list[i], RRLWR_PKE_LOGQ, seedA, RRLWR_PKE_SEED_A_LEN, 0);
    poly_uniform(&b_list[i], RRLWR_PKE_LOGQ, seedB, RRLWR_PKE_SEED_A_LEN, 0);
  }
}

static void bench_poly_ntt(void) {
  poly a, b, r;

  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    a = a_list[i];
    b = b_list[i];

    // poly_ntt32(&a, RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV, rrlwr_pke_zetas);
    poly_ntt32(&b, RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV, rrlwr_pke_zetas);
    poly_basemul32(&r, &a, &b, RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV);
    poly_invntt32(&r, RRLWR_PKE_PRIME, RRLWR_PKE_PRIMEINV,
                  RRLWR_NTTINV_FINALCONST, rrlwr_pke_zetas);
    poly_conditional_final_reduce32(&r, RRLWR_PKE_PRIME);

    t[i] = cpucycles();
  }

  print_results("poly NTT + pointwise + intt/finalreduce", t, NUMBER_OF_TESTS);
}

static void bench_poly_q13_toom(void) {
  poly r;

  for(unsigned int i = 0; i < NUMBER_OF_TESTS; i++) {
    poly_mul_q13(&r, &a_list[i], &b_list[i]);
    t[i] = cpucycles();
  }

  print_results("poly q13 Toom-Cook + negacyclic fold", t, NUMBER_OF_TESTS);
}

int main(void) {
  init_inputs();

  bench_poly_ntt();
  bench_poly_q13_toom();

  return 0;
}
