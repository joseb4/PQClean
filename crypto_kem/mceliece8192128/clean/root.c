/*
  This file is for evaluating a polynomial at one or more field elements
*/

#include "root.h"
#include "gf.h"
#include "params.h"

#include <stdio.h>

/* input: polynomial f and field element a */
/* return f(a) */
<<<<<<< HEAD
gf eval(gf *f, gf a) {
=======
gf eval(gf *f, gf a, const int degree) {
>>>>>>> master
    int i;
    gf r;

    r = f[ SYS_T ];

<<<<<<< HEAD
    for (i = SYS_T - 1; i >= 0; i--) {
=======
    for (i = SYS_T - 1; i >= degree; i--) {
>>>>>>> master
        r = gf_mul(r, a);
        r = gf_add(r, f[i]);
    }

    return r;
}

/* input: polynomial f and list of field elements L */
/* output: out = [ f(a) for a in L ] */
void root(gf *out, gf *f, gf *L) {
    int i;
<<<<<<< HEAD

    for (i = 0; i < SYS_N; i++) {
        out[i] = eval(f, L[i]);
=======
    int degree = -1;
    for (i = 0; i <= SYS_T; i++) {
        if (f[i] != 0) {
            degree = i;
            break;
        }
    }

    for (i = 0; i < SYS_N; i++) {
        out[i] = eval(f, L[i], degree);
>>>>>>> master
    }
}
