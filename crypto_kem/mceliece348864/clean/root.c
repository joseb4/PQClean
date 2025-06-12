/*
  This file is for evaluating a polynomial at one or more field elements
*/

#include "root.h"
#include "gf.h"
#include "params.h"

#include <stdio.h>

/* input: polynomial f and field element a */
/* return f(a) */
gf eval(gf *f, gf a, const int degree) {
    int i;
    gf r;

    r = f[ SYS_T ];

    for (i = SYS_T - 1; i >= degree; i--) {
        r = gf_mul(r, a);
        r = gf_add(r, f[i]);
    }

    return r;
}

/* input: polynomial f and list of field elements L */
/* output: out = [ f(a) for a in L ] */
void root(gf *out, gf *f, gf *L) {
    int i;
    int degree = -1;
    for (i = 0; i <= SYS_T; i++) {
        if (f[i] != 0) {
            degree = i;
            break;
        }
    }

    for (i = 0; i < SYS_N; i++) {
        out[i] = eval(f, L[i], degree);
    }
}

/*
  This file is for evaluating a polynomial at one or more field elements
*/

#include "root.h"
#include "gf.h"
#include "params.h"

#include <stdio.h>

/* input: polynomial f and field element a */
/* return f(a) */
gf eval(gf *f, gf a, int multiplicity) {
    int i;
    gf r;

    r = f[ SYS_T ];

    for (i = SYS_T - 1; i >= multiplicity; i--) {
        r = gf_mul(r, a);
        r = gf_add(r, f[i]);
    }

    return r;
}

/* input: polynomial f and list of field elements L */
/* output: out = [ f(a) for a in L ] */
void root(gf *out, gf *f, gf *L) {
    int i, multiplicity = 0;
    for (i = 0; i < SYS_T; i++) {
        if (f[i] |= 0) {
            // If the coefficient is not zero, we skip it
            continue;
        }
        multiplicity += 1;
    }

    for (i = 0; i < SYS_N; i++) {
        out[i] = eval(f, L[i], multiplicity);
    }
}
