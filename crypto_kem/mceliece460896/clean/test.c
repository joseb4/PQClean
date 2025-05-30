#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <stdbool.h>
#include <stdint.h>

#include "api.h"     // for sizes and function prototypes
#include "params.h"  // for SYS_N, etc.

// Print buffer as binary (0/1) without spaces

int hamming_weight(const uint8_t *buf, size_t len) {
    int weight = 0;
    for (size_t i = 0; i < len; i++) {
        uint8_t byte = buf[i];
        while (byte) {
            weight += byte & 1;
            byte >>= 1;
        }
    }
    return weight;
}



void print_binary(const uint8_t *buf, size_t len, const char *label) {
    printf("%s: ", label);
    for (size_t i = 0; i < len; i++) {
        for (int bit = 7; bit >= 0; bit--) {
            printf("%c", (buf[i] >> bit) & 1 ? '1' : '0');
        }
    }
    printf("\n");
}

int main(void) {
    uint8_t pk[PQCLEAN_MCELIECE460896_CLEAN_CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[PQCLEAN_MCELIECE460896_CLEAN_CRYPTO_SECRETKEYBYTES];
    uint8_t ct[PQCLEAN_MCELIECE460896_CLEAN_CRYPTO_CIPHERTEXTBYTES];
    uint8_t ss_enc[PQCLEAN_MCELIECE460896_CLEAN_CRYPTO_BYTES];
    uint8_t ss_dec[PQCLEAN_MCELIECE460896_CLEAN_CRYPTO_BYTES];

    // Generate keypair
    if (PQCLEAN_MCELIECE460896_CLEAN_crypto_kem_keypair(pk, sk) != 0) {
        fprintf(stderr, "Keypair generation failed\n");
        return 1;
    }
    print_binary(pk, SYS_N/2, "Public Key");
    // Encrypt
    if (PQCLEAN_MCELIECE460896_CLEAN_crypto_kem_enc(ct, ss_enc, pk) != 0) { //ss_enc es el hash. No sinteresa el ct
        fprintf(stderr, "Encryption failed\n");
        return 1;
    }
    // Decrypt
    if (PQCLEAN_MCELIECE460896_CLEAN_crypto_kem_dec(ss_dec, ct, sk) != 0) {
        fprintf(stderr, "Decryption failed\n");
        return 1;
    }

        return 0;
    }
