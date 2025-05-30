#ifndef PK_GEN_H
#define PK_GEN_H
/*
  This file is for public-key generation
*/

#include "namespace.h"

#define pk_gen CRYPTO_NAMESPACE(pk_gen)
#define compute_generator_matrix CRYPTO_NAMESPACE(compute_generator_matrix)
#define gen_G_matrix CRYPTO_NAMESPACE(gen_G_matrix)
#define build_G_matrix CRYPTO_NAMESPACE(build_G_matrix)

#include "gf.h"

int pk_gen(unsigned char *pk, unsigned char *sk, const uint32_t *perm, int16_t *pi);
void compute_generator_matrix(unsigned char *G, const unsigned char *pk);
void gen_G_matrix(uint8_t *G, const uint8_t *pk);
void build_G_matrix(unsigned char G[][SYS_N / 8], const unsigned char *pk);

#endif
