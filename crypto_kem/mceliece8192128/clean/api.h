#ifndef PQCLEAN_MCELIECE8192128_CLEAN_API_H
#define PQCLEAN_MCELIECE8192128_CLEAN_API_H

#include <stdint.h>

#define PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_ALGNAME "Classic McEliece 8192128"
#define PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_PUBLICKEYBYTES PK_NROWS*PK_NCOLS/8
#define PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_SECRETKEYBYTES (32+8)+SYS_T*2+((1 << (GFBITS - 4)) * (2 * GFBITS - 1)) + SYS_N / 8 //6492
#define PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_CIPHERTEXTBYTES PK_NROWS/8 //96 size_of( H*e) = PK_NCOLS
#define PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_BYTES 32

int PQCLEAN_MCELIECE8192128_CLEAN_crypto_kem_enc(
    uint8_t *c,
    uint8_t *key,
    const uint8_t *pk
);

int PQCLEAN_MCELIECE8192128_CLEAN_crypto_kem_dec(
    uint8_t *key,
    const uint8_t *c,
    const uint8_t *sk
);

int PQCLEAN_MCELIECE8192128_CLEAN_crypto_kem_keypair
(
    uint8_t *pk,
    uint8_t *sk
);
int PQCLEAN_MCELIECE8192128_CLEAN_extract_preimage(
    unsigned char *preimage_out, 
    const unsigned char *c, 
    const unsigned char *sk, 
    unsigned char *e
);
void PQCLEAN_MCELIECE8192128_CLEAN_encrypt(
    unsigned char *s, 
    const unsigned char *pk, 
    unsigned char *e
);
void PQCLEAN_MCELIECE8192128_CLEAN_syndrome(
    unsigned char *s, 
    const unsigned char *pk, 
    const unsigned char *e
);
void PQCLEAN_MCELIECE8192128_CLEAN_codeword(
    unsigned char *xG, 
    const unsigned char *pk, 
    const unsigned char *x
);
void PQCLEAN_MCELIECE8192128_CLEAN_gen_weight(
    unsigned char *e, 
    int weight
);

#endif
