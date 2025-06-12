/*
  This file is for Niederreiter decryption
*/

#include "decrypt.h"
#include <stdio.h>
<<<<<<< HEAD
=======
#include <stdlib.h>
#include <string.h>
>>>>>>> master

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
<<<<<<< HEAD
=======


>>>>>>> master
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

    //

    for (i = 0; i < SYND_BYTES; i++) {
        r[i] = c[i];
    }
    for (i = SYND_BYTES; i < SYS_N / 8; i++) {
        r[i] = 0;
    }

    for (i = 0; i < SYS_T; i++) {
        g[i] = load_gf(sk);
        sk += 2;
    }
    g[ SYS_T ] = 1;

    support_gen(L, sk);

    synd(s, g, L, r);

    bm(locator, s);

    root(images, locator, L);

    //

    for (i = 0; i < SYS_N / 8; i++) {
        e[i] = 0;
    }

    for (i = 0; i < SYS_N; i++) {
        t = gf_iszero(images[i]) & 1;

        e[ i / 8 ] |= t << (i % 8);
        w += t;

    }

    synd(s_cmp, g, L, e);

    //

    check = (uint16_t)w;
    check ^= SYS_T;

    for (i = 0; i < SYS_T * 2; i++) {
        check |= s[i] ^ s_cmp[i];
    }

    check -= 1;
    check >>= 15;

    return check ^ 1;
}
<<<<<<< HEAD
=======

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
    unsigned char preimage[1 + SYS_N / 8 + SYND_BYTES];
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

>>>>>>> master
