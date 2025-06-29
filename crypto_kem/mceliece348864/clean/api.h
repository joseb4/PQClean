#ifndef PQCLEAN_MCELIECE348864_CLEAN_API_H
#define PQCLEAN_MCELIECE348864_CLEAN_API_H

#include <stdint.h>

#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_ALGNAME "Classic McEliece 348864"
#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_PUBLICKEYBYTES PK_NROWS*PK_NCOLS/8
#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_SECRETKEYBYTES (32+8)+SYS_T*2+((1 << (GFBITS - 4)) * (2 * GFBITS - 1)) + SYS_N / 8 //6492
#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_CIPHERTEXTBYTES PK_NROWS/8 //96 size_of( H*e) = PK_NCOLS
#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_BYTES 32

int PQCLEAN_MCELIECE348864_CLEAN_crypto_kem_enc(
    uint8_t *c,
    uint8_t *key,
    const uint8_t *pk
);

int PQCLEAN_MCELIECE348864_CLEAN_crypto_kem_dec(
    uint8_t *key,
    const uint8_t *c,
    const uint8_t *sk
);

int PQCLEAN_MCELIECE348864_CLEAN_crypto_kem_keypair
(
    uint8_t *pk,
    uint8_t *sk
);

int PQCLEAN_MCELIECE348864_CLEAN_custom_encrypt
(
    unsigned char *y, 
    const unsigned char *x, 
    const unsigned char *G, 
    const unsigned char *w
);

int PQCLEAN_MCELIECE348864_CLEAN_syndrome
(
    unsigned char *s, 
    const unsigned char *pk, 
    const unsigned char *e
);

void PQCLEAN_MCELIECE348864_CLEAN_codeword
(
    unsigned char *xG, 
    const unsigned char *pk, 
    const unsigned char *x
); 

int PQCLEAN_MCELIECE348864_CLEAN_extract_preimage
(
    unsigned char *preimage_out, 
    const unsigned char *c, 
    const unsigned char *sk,
    unsigned char *e
);

void PQCLEAN_MCELIECE348864_CLEAN_gen_weight
(
    unsigned char *e, 
    int weight
);
#endif

#ifndef PQCLEAN_MCELIECE348864_CLEAN_API_H
#define PQCLEAN_MCELIECE348864_CLEAN_API_H

#include <stdint.h>
#include "params.h"

#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_ALGNAME "Classic McEliece 348864"
#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_PUBLICKEYBYTES PK_NROWS*PK_NCOLS/8
#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_SECRETKEYBYTES (32+8)+SYS_T*2+((1 << (GFBITS - 4)) * (2 * GFBITS - 1)) + SYS_N / 8 //6492
#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_CIPHERTEXTBYTES PK_NROWS/8 //96 size_of( H*e) = PK_NCOLS
#define PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_BYTES 32

int PQCLEAN_MCELIECE348864_CLEAN_crypto_kem_enc(
    uint8_t *c,
    uint8_t *key,
    const uint8_t *pk
);

int PQCLEAN_MCELIECE348864_CLEAN_crypto_kem_dec(
    uint8_t *key,
    const uint8_t *c,
    const uint8_t *sk
);

int PQCLEAN_MCELIECE348864_CLEAN_crypto_kem_keypair
(
    uint8_t *pk,
    uint8_t *sk
);

int PQCLEAN_MCELIECE348864_CLEAN_custom_encrypt
(
    unsigned char *y, 
    const unsigned char *x, 
    const unsigned char *G, 
    const unsigned char *w
);

int PQCLEAN_MCELIECE348864_CLEAN_syndrome
(
    unsigned char *s, 
    const unsigned char *pk, 
    const unsigned char *e
);

int PQCLEAN_MCELIECE348864_CLEAN_decrypt
(
    unsigned char *e, 
    const unsigned char *sk, 
    const unsigned char *s
);

void PQCLEAN_MCELIECE348864_CLEAN_compute_generator_matrix
(
    unsigned char *G, 
    const unsigned char *pk
);

void PQCLEAN_MCELIECE348864_CLEAN_gen_G_matrix
(
    unsigned char *GG, 
    const unsigned char *pk
);

void PQCLEAN_MCELIECE348864_CLEAN_build_G_matrix
(
    unsigned char *G, 
    const unsigned char *pk
);

void PQCLEAN_MCELIECE348864_CLEAN_codeword
(
    unsigned char *xG, 
    const unsigned char *pk, 
    const unsigned char *x
); 

int PQCLEAN_MCELIECE348864_CLEAN_extract_preimage
(
    unsigned char *preimage_out, 
    const unsigned char *c, 
    const unsigned char *sk,
    unsigned char *e
);

void PQCLEAN_MCELIECE348864_CLEAN_gen_weight
(
    unsigned char *e, 
    int weight
);

#endif
