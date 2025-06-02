/*
  This file is for Niederreiter decryption
*/

#include "decrypt.h"
#include <stdio.h>
#include <string.h>
#include "benes.h"
#include "bm.h"
#include "gf.h"
#include "params.h"
#include "root.h"
#include "synd.h"
#include "util.h"

/* Niederreiter decryption with the Berlekamp decoder */
/* intput: sk, secret key */
/*         c, ciphertext */
/* output: e, error vector */
/* return: 0 for success; 1 for failure */
int decrypt(unsigned char *e, const unsigned char *sk, const unsigned char *c) {
    int i, w = 0;
    uint16_t check;

    unsigned char r[ SYS_N / 8 ];

    gf g[ SYS_T + 1 ];
    gf L[ SYS_N ];

    gf s[ SYS_T * 2 ];
    gf s_cmp[ SYS_T * 2 ];
    gf locator[ SYS_T + 1 ];
    gf images[ SYS_N ];

    gf t;

    // C = He
    // extend v = (C,0,...,0) 

    for (i = 0; i < SYND_BYTES; i++) {
        r[i] = c[i];
    }
    for (i = SYND_BYTES; i < SYS_N / 8; i++) {
        r[i] = 0;
    }
    // Find the unique c\in F^n/2 s.t,. 
    //  1) Hc = 0
    //  2) wt(c-v) <= SYS_T
    //  3) e = v + c

    // 1) Load the Goppa polynomial and support
    for (i = 0; i < SYS_T; i++) {
        g[i] = load_gf(sk);
        sk += 2;
    }
    g[ SYS_T ] = 1;
    // Load the alpha
    support_gen(L, sk);
    // Calcualte syndrome of r, s = Syn(r) = 1/z + 1/(1+z) mod g(z). This will be an element in GF(2)[z]/g(z)
    synd(s, g, L, r);
    // Calculate the error locator polynomial with Berlecamp
    bm(locator, s);
    // The roots alpha_i indicate the error
    root(images, locator, L);// When deg\neq SYS_T, there is one extra root because the polynomial looks like x^m*(poly)
    //

    for (i = 0; i < SYS_N / 8; i++) {
        e[i] = 0;
    }

    for (i = 0; i < SYS_N; i++) {
        //printf("i = %d, image = %d\n", i, images[i]);
        t = gf_iszero(images[i]) & 1;
        if (gf_iszero(images[i]) & 1)
            printf("Error: image[%d] is zero\n", i);

        e[ i / 8 ] |= t << (i % 8);
        //e[ (SYS_N -i) / 8 ] |= t << (i % 8);
        
        w += t;

    }
    //for (i = 0; i < 500; i++) 
        //printf("i = %d, image = %d\n", i, images[i]);
    printf("\n%d errors located\n", w);

    synd(s_cmp, g, L, e);

    // This focres the hammight // weight of the error vector to be SYS_T

    if (w > SYS_T) {
        // If the weight is greater than SYS_T, we return an error
        printf("Too many errors located");
        return 1;
    }
    check = (uint16_t)0;
    //check ^= SYS_T;
    // Verify the syndrome
    for (i = 0; i < SYS_T * 2; i++) {
        check |= s[i] ^ s_cmp[i];
    }

    check -= 1;
    check >>= 15;

    return check ^ 1;
}

int extract_preimage(
    unsigned char *pre,
    const unsigned char *c,
    const unsigned char *sk,
    unsigned char *e
) {
    int i;

    unsigned char ret_decrypt = 0;

    uint16_t m;

    //unsigned char e[ SYS_N / 8 ];
    unsigned char preimage[ 1 + SYS_N / 8 + SYND_BYTES ];
    unsigned char *x = preimage;
    const unsigned char *s = sk + 40 + IRR_BYTES + COND_BYTES;

    //

    ret_decrypt = (unsigned char)decrypt(e, sk + 40, c); //Given the cyphertext c, decrypts it using the secret key sk and stores the error vector in e

    m = ret_decrypt;
    m -= 1;
    m >>= 8;

    *x++ = m & 1;
    for (i = 0; i < SYS_N / 8; i++) {
        *x++ = (~m & s[i]) | (m & e[i]);
    }

    for (i = 0; i < SYND_BYTES; i++) {
        *x++ = c[i];
    }
    memcpy(pre, preimage + 1, SYS_N / 8);
    return 0;
}