#ifndef ENCRYPT_H
#define ENCRYPT_H
/*
  This file is for Niederreiter encryption
*/

#include "namespace.h"

#define encrypt CRYPTO_NAMESPACE(encrypt)
#define syndrome CRYPTO_NAMESPACE(syndrome)
#define codeword CRYPTO_NAMESPACE(codeword)
#define gen_weight CRYPTO_NAMESPACE(gen_weight)

void encrypt(unsigned char *s, const unsigned char *pk, unsigned char *e);
void syndrome(unsigned char *s, const unsigned char *pk, const unsigned char *e);
void codeword(unsigned char *xG, const unsigned char *pk, const unsigned char *x);
void gen_weight(unsigned char *e, int weight);


#endif
