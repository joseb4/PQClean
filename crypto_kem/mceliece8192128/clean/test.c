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

    printf("%d\n",PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_PUBLICKEYBYTES);
    printf("%d\n",PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_SECRETKEYBYTES);
    printf("%d\n",PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_CIPHERTEXTBYTES);
    
    printf("KEM Example\n");
    printf("============\n");
    uint8_t pk[PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_SECRETKEYBYTES];
    uint8_t ct[PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_CIPHERTEXTBYTES];
    uint8_t ss_enc[PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_BYTES];
    uint8_t ss_dec[PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_BYTES];

    // Generate keypair
    if (PQCLEAN_MCELIECE8192128_CLEAN_crypto_kem_keypair(pk, sk) != 0) {
        fprintf(stderr, "Keypair generation failed\n");
        return 1;
    }
    
    // Encrypt
    if (PQCLEAN_MCELIECE8192128_CLEAN_crypto_kem_enc(ct, ss_enc, pk) != 0) { //ss_enc es el hash. No sinteresa el ct
        fprintf(stderr, "Encryption failed\n");
        return 1;
    }
    /*
    // Decrypt
    if (PQCLEAN_MCELIECE8192128_CLEAN_crypto_kem_dec(ss_dec, ct, sk) != 0) {
        fprintf(stderr, "Decryption failed\n");
        return 1;
    }
*/
    
    // Print original shared secret key
    print_binary(ss_enc, PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_BYTES, "\nShared Secret (Enc)");

    // Print ciphertext
    print_binary(ct, PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_CIPHERTEXTBYTES, "Ciphertext");

    // Print decrypted shared secret key
    print_binary(ss_dec, PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_BYTES, "Shared Secret (Dec)");

    // Verify keys match
    if (memcmp(ss_enc, ss_dec, PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_BYTES) == 0) {
        printf("✅ Shared secrets match: decryption successful.");
    } else {
        printf("❌ Shared secrets differ: decryption failed.");
    }
    
    
    unsigned char y[SYS_N/8];
    unsigned char x[PK_ROW_BYTES]; // PK_NCOLS bits, packed into bytes
    unsigned char syn[SYND_BYTES];
    unsigned char xG[SYS_N/8];
    unsigned char w[SYS_N/8];
    unsigned char e[SYS_N/8];

    printf("\n\nSECRET ENCRYPTION WITH BIOMETRY y = xG + w. DEMONSTRATOR\n");
    printf("============\n");

    // Generate random message x of correct length and calculate corresponding codeword x*G
    for (size_t i = 0; i < PK_ROW_BYTES; i++) {
            x[i] = rand() % 256;
    }
    print_binary(x, PK_ROW_BYTES, "x (random secret to be encrypted)");
    PQCLEAN_MCELIECE8192128_CLEAN_codeword(xG, pk, x); 
    print_binary(xG, SYS_N/8, "\nCodeword");
    PQCLEAN_MCELIECE8192128_CLEAN_syndrome(syn, pk,xG);
    // Generate encription biometry as random array of bits
    int weight = 63;
    PQCLEAN_MCELIECE8192128_CLEAN_gen_weight(w, weight);
    print_binary(w, SYS_N/8, "\nBiometry error vector w");
    printf("Hamming weight of w: %d\n\n", hamming_weight(w, SYS_N/8));
    PQCLEAN_MCELIECE8192128_CLEAN_syndrome(syn, pk, w);
    print_binary(syn, SYND_BYTES, "Syndrome of error vector w");
    printf("Hammight weight of the syndrome: %d\n\n", hamming_weight(syn, SYND_BYTES));
    
    unsigned char decrypted_w[SYS_N/8];
    PQCLEAN_MCELIECE8192128_CLEAN_extract_preimage(decrypted_w, syn, sk, e);
    
    print_binary(e, SYS_N/8, "\n\nDecrypted error vector e");
    print_binary(decrypted_w, SYS_N/8, "\n\nDecrypted error vector w");
    printf("Hamming weight of decrypted_w: %d\n\n", hamming_weight(decrypted_w, SYS_N/8));
    if (memcmp(e, w, SYS_N/8) == 0) {
        printf("\n✅ Decryption successful: decrypted error vector matches original.\n");
    } else {
        printf("\n❌ Decryption failed: decrypted error vector does not match original.\n");
    }

    PQCLEAN_MCELIECE8192128_CLEAN_syndrome(syn, pk, decrypted_w);
    print_binary(syn, SYND_BYTES, "Syndrome of decrypted error vector w");
        return 0;
    }
