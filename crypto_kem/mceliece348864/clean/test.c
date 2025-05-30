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
    uint8_t pk[PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_SECRETKEYBYTES];
    uint8_t ct[PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_CIPHERTEXTBYTES];
    uint8_t ss_enc[PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_BYTES];
    uint8_t ss_dec[PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_BYTES];

    // Generate keypair
    if (PQCLEAN_MCELIECE348864_CLEAN_crypto_kem_keypair(pk, sk) != 0) {
        fprintf(stderr, "Keypair generation failed\n");
        return 1;
    }
    // Encrypt
    if (PQCLEAN_MCELIECE348864_CLEAN_crypto_kem_enc(ct, ss_enc, pk) != 0) { //ss_enc es el hash. No sinteresa el ct
        fprintf(stderr, "Encryption failed\n");
        return 1;
    }
    // Decrypt
    if (PQCLEAN_MCELIECE348864_CLEAN_crypto_kem_dec(ss_dec, ct, sk) != 0) {
        fprintf(stderr, "Decryption failed\n");
        return 1;
    }


    // Print original shared secret key
    print_binary(ss_enc, PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_BYTES, "Shared Secret (Enc)");

    // Print ciphertext
    print_binary(ct, PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_CIPHERTEXTBYTES, "Ciphertext");

    // Print decrypted shared secret key
    print_binary(ss_dec, PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_BYTES, "Shared Secret (Dec)");

    // Verify keys match
    if (memcmp(ss_enc, ss_dec, PQCLEAN_MCELIECE348864_CLEAN_CRYPTO_BYTES) == 0) {
        printf("✅ Shared secrets match: decryption successful.\n\n");
    } else {
        printf("❌ Shared secrets differ: decryption failed.\n\n");
    }
    
    unsigned char y[SYS_N/8];
    unsigned char x[(PK_NCOLS + 7) / 8]; // PK_NCOLS bits, packed into bytes
    unsigned char G[PK_NROWS * SYS_N];
    unsigned char GG[PK_NCOLS * SYS_N/8];
    unsigned char syn[SYND_BYTES];
    unsigned char xG[SYS_N/8];

    unsigned char ww[SYS_N/8] = {0};
    bool used2[SYS_N] = {0};
    int weight2 = 0;
    while (weight2 < SYS_T) {
        int bitpos2 = rand() % SYS_N;
        if (!used2[bitpos2]) {
            used2[bitpos2] = true;
            ww[bitpos2 / 8] |= 1 << (bitpos2 % 8);
            weight2++;
        }
    }

    // Generate random message x of correct length
    for (size_t i = 0; i < PK_NCOLS; i++) {
        if (rand() & 1) {
            x[i/8] |= 1 << (i % 8);
        }
    }
    print_binary(x, PK_NROWS/8, "x (random vector)");
    PQCLEAN_MCELIECE348864_CLEAN_codeword(xG, pk, x); 
    print_binary(xG, SYS_N/8, "Codeword");
    //PQCLEAN_MCELIECE348864_CLEAN_build_G_matrix(G, pk);
    //print_binary(pk, SYS_N, "H");
    //print_binary(G, SYS_N, "G");
    // Encrypt using custom_encrypt: y = x·G ⊕ ww
    //PQCLEAN_MCELIECE348864_CLEAN_custom_encrypt(y, x, G, ww);
    // Dencription: calculate e= yH^T, from e extract, xG. Then ww = xG ⊕ y
    if (PQCLEAN_MCELIECE348864_CLEAN_syndrome(syn, pk,xG) != 0) {
    fprintf(stderr, "Encryption failed\n");
    return 1;
    }    
    print_binary(syn, SYS_N/8, "syndrome of worderod");
/*


    
    // Decrypt
    PQCLEAN_MCELIECE348864_CLEAN_decrypt(xG, sk, ct);

    uint8_t ww_prime[SYS_N / 8];
    for (size_t i = 0; i < SYS_N / 8; i++) {
        ww_prime[i] = y[i] ^ xG[i];
    }
    // Print original shared secret key
    print_binary(G, PK_NROWS, "G");
    print_binary(x, PK_NROWS, "x");
    print_binary(xG, SYS_N/8, "xG");
    print_binary(ww, SYS_N/8, "Noise");

    // Print ciphertext
    print_binary(ww_prime, SYS_N/8, "Detected noise ");
    
    
    if (memcmp(ww, ww_prime, SYS_N / 8) == 0) {
        printf("✅ Success: w == w'\n");
    } else {
        printf("❌ Error: w != w'\n");
    }

    uint8_t diff[SYS_N / 8];
for (size_t i = 0; i < SYS_N / 8; i++) {
    diff[i] = ww[i] ^ ww_prime[i];
}

// Print Hamming weights
int hw_ww = hamming_weight(ww, SYS_N / 8);
int hw_ww_prime = hamming_weight(ww_prime, SYS_N / 8);
int hw_diff = hamming_weight(diff, SYS_N / 8);

printf("Hamming weight of ww:        %d\n", hw_ww);
printf("Hamming weight of ww_prime:  %d\n", hw_ww_prime);
printf("Hamming weight of difference: %d\n", hw_diff);
        */
        return 0;
    }
