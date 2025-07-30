#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include "api.h"
#include "crypto_hash.h"
#include "params.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



uint8_t patterns[5] = {0x0, 0x1, 0x3, 0x7, 0xF};
const int pattern_indices[16] = {
    [0x0] = 0, [0x1] = 1, [0x3] = 2, [0x7] = 3, [0xF] = 4
};


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
// Function to modify w_bio with a fixed number of changes
void modify_w_bio(unsigned char *w_bio_copy, const unsigned char *w_bio, int target_changes) {
    uint8_t w_bio_diff[SYS_N/8];
    memcpy(w_bio_copy, w_bio, SYS_N/8);
    int total_changes = 0;
    while (total_changes < target_changes) {
        // Randomly select a nibble
        int nibble_idx = rand() % 128;
        int byte_idx = nibble_idx / 2;
        int is_high_nibble = (nibble_idx % 2) == 0;

        // Extract the current nibble
        uint8_t current_nibble;
        if (is_high_nibble) {
            current_nibble = (w_bio_copy[byte_idx] >> 4) & 0xF;
        } else {
            current_nibble = w_bio_copy[byte_idx] & 0xF;
        }

        // Find current pattern index (default to 0 if invalid)
        int current_idx = 0;
        for (int i = 0; i < 5; i++) {
            if (patterns[i] == current_nibble) {
                current_idx = i;
                break;
            }
        }
        int direction = (rand() % 2) ? 1 : -1;  // Randomly choose +1 or -1
        if (current_idx+direction < 0 || current_idx+direction >= 5) {
            // If out of bounds, reverse direction
            direction *= -1;
        }
        int new_idx = current_idx + direction;
        // Update the nibble in w_bio
        uint8_t new_nibble = patterns[new_idx];
        if (is_high_nibble) {
            w_bio_copy[byte_idx] = (w_bio_copy[byte_idx] & 0x0F) | (new_nibble << 4);
        } else {
            w_bio_copy[byte_idx] = (w_bio_copy[byte_idx] & 0xF0) | new_nibble;
        }


        for (int i = 0; i < SYS_N/8; i++)
            w_bio_diff[i] = w_bio_copy[i] ^ w_bio[i];
        total_changes = hamming_weight(w_bio_diff, SYS_N/8);
    }
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

void calc_wpass(unsigned char *w_pass, char *password, uint8_t *id_x) {
    size_t seed_len = strlen(password) + 32; // id_x is 32 bytes (SHA-256)
    uint8_t *seed = malloc(seed_len);
    memcpy(seed, password, strlen(password));
    memcpy(seed + strlen(password), id_x, 32);
    shake(w_pass, SYS_N/8, seed, seed_len);
    free(seed);
}

int main(void) {
    srand(time(NULL)); // Seed randomness

    uint8_t pk[PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[PQCLEAN_MCELIECE8192128_CLEAN_CRYPTO_SECRETKEYBYTES];
    unsigned char v[SYS_N/8], x[PK_ROW_BYTES], syn[SYND_BYTES];
    unsigned char xG[SYS_N/8], w_bio[SYS_N/8] = {0}, w_bio_copy[SYS_N/8], w_pass[SYS_N/8], w_pass_copy[SYS_N/8], w[SYS_N/8];
    uint8_t id_x[32]; // SHA-256 hash of x

    // Generate keypair
    if (PQCLEAN_MCELIECE8192128_CLEAN_crypto_kem_keypair(pk, sk) != 0) {
        fprintf(stderr, "Keypair generation failed\n");
        return 1;
    }

    // Generate w_bio (fixed for all iterations)
    for (size_t i = 0; i < 128; i++) {
        uint8_t chosen = patterns[rand() % 5];
        size_t byte_idx = i / 2;
        if (i % 2 == 0) {
            w_bio[byte_idx] |= (chosen << 4); // high nibble
        } else {
            w_bio[byte_idx] |= chosen;        // low nibble
        }
    }

    char *password = "example_password";

    for (int iter = 0; iter < 10000; iter++) {
        // Generate random x
        for (size_t i = 0; i < PK_ROW_BYTES; i++) {
            x[i] = rand() % 256;
        }
        shake(id_x, 32, x, PK_ROW_BYTES); // Hash x to id_x
        calc_wpass(w_pass, password, id_x);

        PQCLEAN_MCELIECE8192128_CLEAN_codeword(xG, pk, x);


        // FOR THE SAME SAMPLE
        for (size_t i = 0; i < SYS_N/8; i++) {
            w[i] = w_bio[i] ^ w_pass[i];
        }
        for (size_t i = 0; i < SYS_N/8; i++) {
            v[i] = w[i] ^ xG[i];
        }
        FILE *f = fopen("same_sample.bin", "ab");
        if (f) {
            fwrite(v, 1, SYS_N/8, f);
            fclose(f);
        } else {
            fprintf(stderr, "Error opening file\n");
        }


        // FOR THE SAME INSTANCE
        double u1 = ((double)rand() + 1) / ((double)RAND_MAX + 2);
        double u2 = ((double)rand() + 1) / ((double)RAND_MAX + 2);
        double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        int changes = (int)(z0 * 10.0 + 70.0);
        if (changes < 0) changes = 0;
        if (changes > 128) changes = 128;
        modify_w_bio(w_bio_copy, w_bio, changes);
        for (size_t i = 0; i < SYS_N/8; i++) {
            w[i] = w_bio_copy[i] ^ w_pass[i];
        }

        for (size_t i = 0; i < SYS_N/8; i++) {
            v[i] = w[i] ^ xG[i];
        }
        FILE *ff = fopen("same_instance.bin", "ab");
        if (ff) {
            fwrite(v, 1, SYS_N/8, ff);
            fclose(ff);
        } else {
            fprintf(stderr, "Error opening file\n");
        }
    

        // FOR RANDOM SAMPLE

    // Generate w_bio (fixed for all iterations)
        for (size_t i = 0; i < 128; i++) {
        uint8_t chosen = patterns[rand() % 5];
        size_t byte_idx = i / 2;
        if (i % 2 == 0) {
            w_bio_copy[byte_idx] |= (chosen << 4); // high nibble
        } else {
            w_bio_copy[byte_idx] |= chosen;        // low nibble
        }
    }
    char *new_pass = "A" + (random() % 26);
    calc_wpass(w_pass_copy, new_pass, id_x);

        for (size_t i = 0; i < SYS_N/8; i++) {
            w[i] = w_bio_copy[i] ^ w_pass[i];
        }

        PQCLEAN_MCELIECE8192128_CLEAN_codeword(xG, pk, x);
        for (size_t i = 0; i < SYS_N/8; i++) {
            v[i] = w[i] ^ xG[i];
        }

        FILE *fff = fopen("diff_sample.bin", "ab");
        if (fff) {
            fwrite(v, 1, SYS_N/8, fff);
            fclose(fff);
        } else {
            fprintf(stderr, "Error opening file\n");
        }
    

    }
    return 0;
}