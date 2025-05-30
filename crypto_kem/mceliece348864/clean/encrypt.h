#ifndef ENCRYPT_H
#define ENCRYPT_H
/*
  This file is for Niederreiter encryption
*/

#include "namespace.h"

#define encrypt CRYPTO_NAMESPACE(encrypt)
#define custom_encrypt CRYPTO_NAMESPACE(custom_encrypt)
#define syndrome CRYPTO_NAMESPACE(syndrome)
#define codeword CRYPTO_NAMESPACE(codeword)

void encrypt(unsigned char *s, const unsigned char *pk, unsigned char *e);
void custom_encrypt(unsigned char *y, const unsigned char *x, const unsigned char *G, const unsigned char *w);
void syndrome(unsigned char *s, const unsigned char *pk, const unsigned char *e);
void codeword(unsigned char *xG, const unsigned char *pk, const unsigned char *x);
#endif
