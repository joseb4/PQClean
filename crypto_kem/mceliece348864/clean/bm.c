/*
  This file is for the Berlekamp-Massey algorithm
  see http://crypto.stanford.edu/~mironov/cs359/massey.pdf
*/

#include "bm.h"
#include "params.h"
#include <stdio.h>

#define min(a, b) (((a) < (b)) ? (a) : (b))

/* the Berlekamp-Massey algorithm */
/* input: s, sequence of field elements */
/* output: out, minimal polynomial of s */
void bm(gf *out, gf *s) {
    int i;

    uint16_t N = 0;
    uint16_t L = 0;
    uint16_t mle;
    uint16_t mne;

    // C: Current minimal polynomial (updated in each step)
    // B: Previous minimal polynomial (updated in each step)
    // T: Temporary storage for the current minimal polynomial
    gf T[ SYS_T + 1  ];
    gf C[ SYS_T + 1 ];
    gf B[ SYS_T + 1 ];

    gf b = 1, d, f;

    // Initialize the polynomials

    for (i = 0; i < SYS_T + 1; i++) {
        C[i] = B[i] = 0;
    }

    // Start with C(x) = 1, B(x) = x
    B[1] = C[0] = 1;

    // Loop for n = 0 to 2T-1
    for (N = 0; N < 2 * SYS_T; N++) {
        // Compute discrepancy d between the current sequence s[x] and the polynomial C(x)
        // if d = 0, C(x) is consistent, continue
        // if d != 0, update the polynomials
        d = 0;
        for (i = 0; i <= min(N, SYS_T); i++) {
            d ^= gf_mul(C[i], s[ N - i]);
        }

        mne = d;
        mne -= 1;
        mne >>= 15;
        mne -= 1;
        // mne = 0xFFFF if d /= d, else mne = 0
        mle = N;
        mle -= 2 * L;
        mle >>= 15;
        mle -= 1;
        mle &= mne;
        // mle = 0xFFFF if N > 2L, else

        // Save polynomial
        for (i = 0; i <= SYS_T; i++) {
            T[i] = C[i];
        }
        // Scalling factor f = d / b
        f = gf_frac(b, d);

        for (i = 0; i <= SYS_T; i++) {
            C[i] ^= gf_mul(f, B[i]) & mne;
        }

        L = (L & ~mle) | ((N + 1 - L) & mle);

        for (i = 0; i <= SYS_T; i++) {
            B[i] = (B[i] & ~mle) | (T[i] & mle);
        }

        b = (b & ~mle) | (d & mle);

        for (i = SYS_T; i >= 1; i--) {
            B[i] = B[i - 1];
        }
        B[0] = 0;
    }
    printf("Coefficients of the error locator polynomial:\n");
    for (i = 0; i <= SYS_T; i++) {
        out[i] = C[ SYS_T - i ];
        //out[i] = C[  i ];
        printf("%d ", out[i]);
    }
}
