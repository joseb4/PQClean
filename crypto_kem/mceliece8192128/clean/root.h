#ifndef ROOT_H
#define ROOT_H
/*
  This file is for evaluating a polynomial at one or more field elements
*/

#include "namespace.h"

#define eval CRYPTO_NAMESPACE(eval)
#define root CRYPTO_NAMESPACE(root)

#include "gf.h"

<<<<<<< HEAD
gf eval(gf *f, gf a);
=======
gf eval(gf *f, gf a, const int degree);
>>>>>>> master
void root(gf *out, gf *f, gf *L);

#endif
