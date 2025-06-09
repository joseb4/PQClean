/*
  This file is for functions for field arithmetic
*/

#include "gf.h"
#include <stddef.h>
#include <stdio.h>
#include "params.h"

#if GFBITS == 13
static const int RED_POLY_EXP[] = {4, 3, 1, 0}; // x^13 + x^4 + x^3 + x + 1
#elif GFBITS == 12
static const int RED_POLY_EXP[] = {3, 0}; // x^12 + x^3 + 1
#else
#error "Unsupported GFBITS value"
#endif
#define RED_POLY_TERMS (sizeof(RED_POLY_EXP)/sizeof(RED_POLY_EXP[0]))

gf gf_iszero(gf a) {
    uint32_t t = a;

    t -= 1;
    t >>= 19;

    return (gf) t;
}

gf gf_add(gf in0, gf in1) {
    return in0 ^ in1;
}

static void reducer(uint64_t *tmp){
    /* The | separates the degree of the polynomial
    0000 0000 0000 000|0 0000 0000 0000
    The final product must be below that line
    Dislapcements are made with the degree of the monomials forming the polynomial, f.i.,
    for f(x) =  x^13+x^4+x^3+x^2+1
    t = tmp & 0x1FF0000;
    tmp ^= (t >> 9) ^ (t >> 10) ^ (t >> 12) ^ (t >> 13);
    t = tmp & 0x000E000;
    tmp ^= (t >> 9) ^ (t >> 10) ^ (t >> 12) ^ (t >> 13);
    since the smallest displacement is 9, the blocksize we pick muset be 9, 1FF or 7FC
    */
    int i;
    size_t j;
    uint64_t t, mask;
    int bits_to_reduce = 64-GFBITS, blocksize = GFBITS-RED_POLY_EXP[0];
    int n_blocks = bits_to_reduce / blocksize;
    int remain = bits_to_reduce - blocksize*n_blocks;


    // Reduce the remaining bits
    mask = ((1ULL<<remain)-1)<<(GFBITS+blocksize*n_blocks); // from the 16 on: 1111 1111 1111 1111 1111 1111 1111 1000 0000 0000 0000
    t = *tmp & mask; 
    for (j = 0; j < RED_POLY_TERMS; j++)
        *tmp ^= (t >> (GFBITS-RED_POLY_EXP[j]));
    // Reduce the rest of the blocks
    for (i = n_blocks-1; i>=0; --i){
        mask = ((1ULL << blocksize)-1)<< (GFBITS+blocksize*i); // from the 13 to 16
        t = *tmp & mask; 
        for (j = 0; j < RED_POLY_TERMS; j++)
            *tmp ^= (t >> (GFBITS-RED_POLY_EXP[j]));
    }

}

gf gf_mul(gf in0, gf in1) {
    int i;
    
    //size_t j;
    //uint64_t mask;
    uint64_t tmp;
    uint64_t t0;
    uint64_t t1;
    //uint64_t t;

    t0 = in0;
    t1 = in1;

    tmp = t0 * (t1 & 1);

    for (i = 1; i < GFBITS; i++) {
        tmp ^= (t0 * (t1 & ((uint64_t)1 << i))); // Multiply coeff by coeff
    }



  
    reducer(&tmp);
    return tmp & GFMASK;
}

/* input: field element in */
/* return: (in^2)^2 */
static inline gf gf_sq2(gf in) {


    const uint64_t B[] = {0x1111111111111111,
                          0x0303030303030303,
                          0x000F000F000F000F,
                          0x000000FF000000FF
                         };
    uint64_t x = in;

    // Double squaring
    x = (x | (x << 24)) & B[3];
    x = (x | (x << 12)) & B[2];
    x = (x | (x << 6)) & B[1];
    x = (x | (x << 3)) & B[0];

    reducer(&x);

    return x & GFMASK;
}

/* input: field element in, m */
/* return: (in^2)*m */
static inline gf gf_sqmul(gf in, gf m) {

    uint64_t x;
    uint64_t t0;
    uint64_t t1;

    t0 = in;
    t1 = m;

    x = (t1 << 6) * (t0 & (1 << 6));
    t0 ^= (t0 << 7);

    x ^= (t1 * (t0 & (0x04001)));
    x ^= (t1 * (t0 & (0x08002))) << 1;
    x ^= (t1 * (t0 & (0x10004))) << 2;
    x ^= (t1 * (t0 & (0x20008))) << 3;
    x ^= (t1 * (t0 & (0x40010))) << 4;
    x ^= (t1 * (t0 & (0x80020))) << 5;

    reducer(&x);

    return x & GFMASK;
}


/* input: field element in, m */
/* return: ((in^2)^2)*m */
static inline gf gf_sq2mul(gf in, gf m) {


    uint64_t x;
    uint64_t t0;
    uint64_t t1;


    t0 = in;
    t1 = m;

    x = (t1 << 18) * (t0 & (1 << 6));

    t0 ^= (t0 << 21);

    x ^= (t1 * (t0 & (0x010000001)));
    x ^= (t1 * (t0 & (0x020000002))) << 3;
    x ^= (t1 * (t0 & (0x040000004))) << 6;
    x ^= (t1 * (t0 & (0x080000008))) << 9;
    x ^= (t1 * (t0 & (0x100000010))) << 12;
    x ^= (t1 * (t0 & (0x200000020))) << 15;

    reducer(&x);
    return x & GFMASK;
}

/* input: field element in */
/* return: in^2 */
static inline gf gf_sq(gf in) {
    const uint64_t B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};

    uint64_t x = in;
    //uint32_t t;

    x = (x | (x << 8)) & B[3];
    x = (x | (x << 4)) & B[2];
    x = (x | (x << 2)) & B[1];
    x = (x | (x << 1)) & B[0];


    reducer(&x);

    return x & ((1 << GFBITS) - 1);
}





/* input: field element den, num */
/* return: (num/den) */
gf gf_frac(gf den, gf num) {
    #if GFBITS == 13
    gf tmp_11;
    gf tmp_1111;
    gf out;

    tmp_11 = gf_sqmul(den, den); // ^11
    tmp_1111 = gf_sq2mul(tmp_11, tmp_11); // ^1111
    out = gf_sq2(tmp_1111);
    out = gf_sq2mul(out, tmp_1111); // ^11111111
    out = gf_sq2(out);
    out = gf_sq2mul(out, tmp_1111); // ^111111111111
    #elif GFBITS == 12
    gf tmp_11;
    gf tmp_1111;

    gf out = den;

    out = gf_sq(out);
    tmp_11 = gf_mul(out, den); // 11

    out = gf_sq(tmp_11);
    out = gf_sq(out);
    tmp_1111 = gf_mul(out, tmp_11); // 1111

    out = gf_sq(tmp_1111);
    out = gf_sq(out);
    out = gf_sq(out);
    out = gf_sq(out);
    out = gf_mul(out, tmp_1111); // 11111111

    out = gf_sq(out);
    out = gf_sq(out);
    out = gf_mul(out, tmp_11); // 1111111111

    out = gf_sq(out);
    out = gf_mul(out, den);
    #endif
    return gf_sqmul(out, num); // ^1111111111110 = ^-1
}

gf gf_inv(gf den) {
    return gf_frac(den, ((gf) 1));
}

/* input: in0, in1 in GF((2^m)^t)*/
/* output: out = in0*in1 */
void GF_mul(gf *out, gf *in0, gf *in1) {
    int i, j;

    gf prod[ SYS_T * 2 - 1 ];

    for (i = 0; i < SYS_T * 2 - 1; i++) {
        prod[i] = 0;
    }

    for (i = 0; i < SYS_T; i++) {
        for (j = 0; j < SYS_T; j++) {
            prod[i + j] ^= gf_mul(in0[i], in1[j]);
        }
    }

    //

    for (i = (SYS_T - 1) * 2; i >= SYS_T; i--) {
        #if GFBITS == 13
        prod[i - SYS_T + 7] ^= prod[i];
        prod[i - SYS_T + 2] ^= prod[i];
        prod[i - SYS_T + 1] ^= prod[i];
        prod[i - SYS_T + 0] ^= prod[i];
        #elif GFBITS == 12
        prod[i - SYS_T + 3] ^= prod[i];
        prod[i - SYS_T + 1] ^= prod[i];
        prod[i - SYS_T + 0] ^= gf_mul(prod[i], (gf) 2);
        #endif
    }

    for (i = 0; i < SYS_T; i++) {
        out[i] = prod[i];
    }
}



