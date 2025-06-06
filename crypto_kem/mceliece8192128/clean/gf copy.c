/*
  This file is for functions for field arithmetic
*/
#include <stddef.h>
#include "gf.h"
#include "params.h"


#if GFBITS == 13
static const int RED_POLY_EXP[] = {13, 4, 3, 1, 0}; // x^13 + x^4 + x^3 + x + 1
#elif GFBITS == 12
static const int RED_POLY_EXP[] = {12, 3, 0}; // x^12 + x^3 + 1
#elif GFBITS == 11
static const int RED_POLY_EXP[] = {11, 2, 0}; // x^12 + x^3 + 1
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

gf gf_mul(gf in0, gf in1) {
    int i;
    


    uint64_t tmp;
    uint64_t t0;
    uint64_t t1;
    uint64_t t;

    t0 = in0;
    t1 = in1;

    tmp = t0 * (t1 & 1);

    for (i = 1; i < GFBITS; i++) {
        tmp ^= (t0 * (t1 & ((uint64_t)1 << i))); // Multiply coeff by coeff
    }

    #if GFBITS == 13
        size_t j;
    uint64_t mask;
    mask = ((1ULL<<(64-GFBITS))-1)<<(GFBITS+2); // from the 16 on: 1111 1111 1111 1111 1111 1111 1111 1000 0000 0000 0000
    

    t = tmp & mask; 
    for (j = 0; j < RED_POLY_TERMS; j++)
        tmp ^= (t >> (GFBITS-RED_POLY_EXP[j]));
    mask = ((1 << 3)-1)<< GFBITS; // from the 13 to 16
    t = tmp & mask; 
    for (j = 0; j < RED_POLY_TERMS; j++)
        tmp ^= (t >> (GFBITS-RED_POLY_EXP[j]));
    return tmp & GFMASK;
    #else
    t = tmp & 0x7FC000;
    tmp ^= t >> 9;
    tmp ^= t >> 12;

    t = tmp & 0x3000;
    tmp ^= t >> 9;
    tmp ^= t >> 12;

    return tmp & ((1 << GFBITS) - 1);
    #endif
}

/* input: field element in */
/* return: in^2 */
gf gf_sq(gf in) {

    const uint64_t B[] = {0x0000000055555555,
                          0x0000000033333333,
                          0x000000000F0F0F0F,
                          0x0000000000FF00FF
                         };
    //const uint32_t B[] = {0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF};

    uint64_t x = in;
    uint64_t t;

    x = (x | (x << 8)) & B[3];
    x = (x | (x << 4)) & B[2];
    x = (x | (x << 2)) & B[1];
    x = (x | (x << 1)) & B[0];

    # if GBITS == 13
    int i;
    size_t j;
    uint64_t mask;
   int bits_to_reduce = 64-GFBITS, chunck = 7;
   mask = ((1ULL << (bits_to_reduce%chunck)) -1)<< ((GFBITS)+chunck*(bits_to_reduce/chunck));
    t = x & mask;
    for (j = 0; j < RED_POLY_TERMS; j++)
        x ^= (t >> (GFBITS-RED_POLY_EXP[j]));

    for (i = (bits_to_reduce/chunck)-1; i >= 0; --i) { //Number of 8 bits chungs
        mask = ((1ULL << chunck) -1)<< ((GFBITS)+chunck*i);
        t = x & mask;
        for (j = 0; j < RED_POLY_TERMS; j++)
            x ^= (t >> (GFBITS-RED_POLY_EXP[j]));
    }
    #else
    t = x & 0x0000007FC000;
    x ^= t >> 9;
    x ^= t >> 12;

    t = x & 0x000000003000;
    x ^= t >> 9;
    x ^= t >> 12;
    #endif


    return x & ((1 << GFBITS) - 1);
}

/* input: field element in */
/* return: (in^2)^2 */
static inline gf gf_sq2(gf in) {
    // In GF(2)(x), if f(x)= a_0 +a_1x + a_2x^2 + ... + a_{m-1}x^{m-1} + a_mx^m, then f(x)^2 = a_0^2 + a_1^2x^2 + a_2^2x^4 + ... + a_{m-1}^2x^{2(m-1)} + a_m^2x^{2m}
    // f(x)^4= (f(x)^2)^2 = a_0^4 + a_1^4x^4 + a_2^4x^8 + ... + a_{m-1}^4x^{4(m-1)} + a_m^4x^{4m}
    int i;
    //uint64_t mask;

    const uint64_t B[] = {0x1111111111111111,
                          0x0303030303030303,
                          0x000F000F000F000F,
                          0x000000FF000000FF
                         };
    // 0001 0001 0001 0001 0001 0001 0001 0001 0001 0001 0001 0001 0001 0001 0001 0001 // coeffs are separated 4 spsaces.
    // 0000 0011 0000 0011 0000 0011 0000 0011 0000 0011 0000 0011 0000 0011 0000 0011
    // 0000 0000 0000 1111 0000 0000 0000 1111 0000 0000 0000 1111 0000 0000 0000 1111
    // 0000 0000 0000 0000 0000 0000 1111 1111 0000 0000 0000 0000 0000 0000 1111 1111  
    //const uint64_t M[] = {0x0001FF0000000000,
    //                      0x000000FF80000000,
    //                      0x000000007FC00000,
    //                      0x00000000003FE000
    //                     };
    // 0000 0000 0000 0001 1111 1111 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000
    // 0000 0000 0000 0000 0000 0000 1111 1111 1000 0000 0000 0000 0000 0000 0000 0000
    // 0000 0000 0000 0000 0000 0000 0000 0000 0111 1111 1100 0000 0000 0000 0000 0000
    // 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0011 1111 1110 0000 0000 0000

    uint64_t x = in;

    uint64_t t;
    size_t j;
    // Gradually sparate the coefficients of the polynomial
    x = (x | (x << 24)) & B[3]; // Displace and take the first 8 bits and the last 8 bits separately
    x = (x | (x << 12)) & B[2];
    x = (x | (x << 6))  & B[1];
    x = (x | (x << 3))  & B[0];
    /*
    ii =37;
    while ( (ii-9) > GFBITS){
        mask = ((1ULL << 9) - 1) << ii;
        ii = ii - 9;
        
        t = x & mask;
        for (j = 0; j < RED_POLY_TERMS; j++)
            x ^= (t >> (GFBITS-RED_POLY_EXP[j]));
    }
    */
   // Reduce the resulting polynomial. 
   // 64 - GFBITS number of bits to displace. We take chuncks of (GFBITS-2), 
   // Reduce remaining bits. There are (64-GFBITS) to reduce in total. We reduce 10 by then, there will be (64-GFBITS)%10 remaining
    uint64_t mask;
   int bits_to_reduce = 64-GFBITS, chunck = 7;
   mask = ((1ULL << (bits_to_reduce%chunck)) -1)<< ((GFBITS)+chunck*(bits_to_reduce/chunck));
    t = x & mask;
    for (j = 0; j < RED_POLY_TERMS; j++)
        x ^= (t >> (GFBITS-RED_POLY_EXP[j]));

    for (i = (bits_to_reduce/chunck)-1; i >= 0; --i) { //Number of 8 bits chungs
        mask = ((1ULL << chunck) -1)<< ((GFBITS)+chunck*i);
        t = x & mask;
        for (j = 0; j < RED_POLY_TERMS; j++)
            x ^= (t >> (GFBITS-RED_POLY_EXP[j]));
    }

    return x & GFMASK;
}

/* input: field element in, m */
/* return: (in^2)*m */
static inline gf gf_sqmul(gf in, gf m) {
    // if in(x) = a_0 + a_1x + a_2x^2 + ... + a_{m-1}x^{m-1} + a_mx^m, then in^2(x) = a_0^2 + a_1^2x^2 + a_2^2x^4 + ... + a_{m-1}^2x^{2(m-1)} + a_m^2x^{2m}
    // then in²(x) * m(x) = (a_0^2 + a_1^2x^2 + a_2^2x^4 + ... + a_{m-1}^2x^{2(m-1)} + a_m^2x^{2m}) * (b_0 + b_1x + b_2x^2 + ... + b_{m-1}x^{m-1} + b_mx^m)
    // = // = a_0^2*b_0 + a_0*b_1x + (a_0^2*b_2 + a_1^2*b_0)x^2 + (a_0^2*b_3 + a_1^2*b_1 + a_2^2*b_0)x^3 + ... + (a_m^2*b_m)x^{2m}
    // = a0*(b0 + b2x^2 + b4x^4 + ... + b2m*x^{2m}) + a1*(b1x + b3x^3 + ... + b(2m-1)x^{2m-1}) + a2*(b2x^2 + b6x^6 + ... + b(2m)x^{2m}) + ... + am^2*(bmx^m)
    int i;
    size_t j;
    uint64_t x;
    uint64_t t0;
    uint64_t t1;
    uint64_t t;

    //const uint64_t M[] = {0x0000001FF0000000, // 0000 0001 1111 1111 0000 0000 0000 0000 0000 
    //                      0x000000000FF80000, // 0000 0000 0000 0000 1111 1111 1000 0000 0000 
    //                      0x000000000007E000  // 0000 0000 0000 0000 0000 0000 0111 1110 0000 
    //                     };

    t0 = in;
    t1 = m;

    x = (t1 << 6) * (t0 & (1 << 6)); //shift 12 if the 6th bit is 1.
    t0 ^= (t0 << 7); // t0 = t0 + t0*x^7, t1*t0 > t1^t0+t1*t0*x^7
    // Shift is due to squaring. For coef 0 no displacement. For a_1 one displacemente, ...
    x ^= (t1 * (t0 & (0x04001)));       // 0000 0100 0000 0000 0001 // Coeffiente 14 was before in 7
    x ^= (t1 * (t0 & (0x08002))) << 1;  // 0000 1000 0000 0000 0010 // multiply coeff 1 to t0, with t1, and siplace 1
    x ^= (t1 * (t0 & (0x10004))) << 2;  // 0001 0000 0000 0000 0100 // multiply coeff 2 to t0, with t1, and siplace 2
    x ^= (t1 * (t0 & (0x20008))) << 3;  // 0010 0000 0000 0000 1000
    x ^= (t1 * (t0 & (0x40010))) << 4;  // 0100 0000 0000 0001 0000
    x ^= (t1 * (t0 & (0x80020))) << 5;  // 1000 0000 0000 0010 0000
    //x = in*in*m;
    uint64_t mask;
   int bits_to_reduce = 64-GFBITS, chunck = 7;
   mask = ((1ULL << (bits_to_reduce%chunck)) -1)<< ((GFBITS)+chunck*(bits_to_reduce/chunck));
    t = x & mask;
    for (j = 0; j < RED_POLY_TERMS; j++)
        x ^= (t >> (GFBITS-RED_POLY_EXP[j]));

    for (i = (bits_to_reduce/chunck)-1; i >= 0; --i) { //Number of 8 bits chungs
        mask = ((1ULL << chunck) -1)<< ((GFBITS)+chunck*i);
        t = x & mask;
        for (j = 0; j < RED_POLY_TERMS; j++)
            x ^= (t >> (GFBITS-RED_POLY_EXP[j]));
    }
    return x & GFMASK;
}

/* input: field element in, m */
/* return: ((in^2)^2)*m */
static inline gf gf_sq2mul(gf in, gf m) {
    int i;
    size_t j;
    uint64_t x;
    uint64_t t0;
    uint64_t t1;
    uint64_t t;

    //const uint64_t M[] = {0x1FF0000000000000,
    //                      0x000FF80000000000,
    //                      0x000007FC00000000,
    //                      0x00000003FE000000,
    //                      0x0000000001FE0000,
    //                      0x000000000001E000 // This is in the 13th position
    //                     };

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

    uint64_t mask;
   int bits_to_reduce = 64-GFBITS, chunck = 7;
   mask = ((1ULL << (bits_to_reduce%chunck)) -1)<< ((GFBITS)+chunck*(bits_to_reduce/chunck));
    t = x & mask;
    for (j = 0; j < RED_POLY_TERMS; j++)
        x ^= (t >> (GFBITS-RED_POLY_EXP[j]));

    for (i = (bits_to_reduce/chunck)-1; i >= 0; --i) { //Number of 8 bits chungs
        mask = ((1ULL << chunck) -1)<< ((GFBITS)+chunck*i);
        t = x & mask;
        for (j = 0; j < RED_POLY_TERMS; j++)
            x ^= (t >> (GFBITS-RED_POLY_EXP[j]));
    }

    return x & GFMASK;
}

/* input: field element den, num */
/* return: (num/den) */

gf gf_frac(gf den, gf num) {
    // den^{-1} = den^{2^q-2} = den^{2^q-1} * den^{-1} = den^{2^q-1} * den^{2^q-2}
    gf tmp_11;
    gf tmp_1111;
    gf out;
    #if GFBITS ==12
        out = den;
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
        out = gf_mul(out, den); // 11111111111 = 2^12-1

        return gf_sqmul(out, num); // den^111111111110*num= num/den
    #else
    tmp_11 = gf_sqmul(den, den); // den^3 = den^2 * den = den^(2^2-1)
    if (GFBITS == 3){ // We want out = den^(2^2-1)
        out = tmp_11; 
        return gf_sqmul(out, num); // returns num * den^{-1}
    }
    if (GFBITS == 4){ //We want out = den^(2^3-1)
        out = gf_sqmul(tmp_11, den); // den^3^2 * den = den^(2^3-2) * den = den^(2^3-1)
        return gf_sqmul(out, num); // returns num * den^{-1}
    }
    tmp_1111 = gf_sq2mul(tmp_11, tmp_11); // den^15 = den^(2^2-1)^(2^2) * den^(2^2-1) = den^(2^4-2^2) * den^(2^2-1) = den^(2^4-1) 
    out = tmp_1111; // den^15 = den^(2^4-1)
    if (GFBITS == 5){
        return gf_sqmul(out, num); // returns num * den^{-1}
    }
    
    if (GFBITS > 5) {
        for (int i = 0; i < (GFBITS-5)/4; i++) { // we want out = den^(2^(GFBITS-1)-1)
            // We already reached it for GFBITS =5. Now we increase by 4 every 4 GFBITS.
            out = gf_sq2(out); // den^60 = (den^(2^4-1))^(2^2) = den^(2^6-2^2)
            out = gf_sq2mul(out, tmp_1111); // den^255 = den^(2^6-2^2)^(2^2) * den^(2^4-1) = den^(2^8-2^4) * den^(2^4-1) = den^(2^8-1)
            // den^(2^(4*floor((GFBITS-1)/4)) - 1)
        }
        // if floor((GFBITS-1)/4)) < (GFBITS-1)/4, We need to add the difference
        if ((GFBITS-1)%4 >1 ){ // If the differnece is 2 or 3, firstly add 2
            // For isntance, if GFBITS = 7,8 den^(2^4-1)^(2^2)*den^(2^2-1)  = den^(2^6-1)
            out = gf_sq2mul(out, tmp_11); // sq2mult to increase by 2, and tmp_11 to ease 2^2
            if ((GFBITS-1)%4 == 2) 
                return gf_sqmul(out, num);
            // If the difference is 3, we need to add one more. F.i. if GBITS = 8, den^(2^6-1)^(2^1)*den^(1) = den^(2^7-1)
            out = gf_sqmul(out, den); // den^(2^1)^2 * den^(2^(4*(GFBITS/4)+2) - 1)= den^(2^(4*(GFBITS/4)+2+1) - 1)
            return gf_sqmul(out, num); // returns num * den^{-1}
        }
        if ((GFBITS-1)%4 == 1){ // For example, if GBITS = 6, den^(2^4-1)^(2^1)*den^(1) = den^(2^5-2+1)  
            out = gf_sqmul(out, den); // Increase by 2  den^(2^(4*(floor(GFBITS-1)/4)) - 1)^(2 ^1)* = den^(2^(4*(GFBITS/4)+1)-2^2 - 1)
            return gf_sqmul(out, num);
        }
    }
    return gf_sqmul(out, num); // returns num * den^{-1}
    #endif
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
