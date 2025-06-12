#ifndef DECRYPT_H
#define DECRYPT_H
/*
  This file is for Nieddereiter decryption
*/

#include "namespace.h"

#define decrypt CRYPTO_NAMESPACE(decrypt)
<<<<<<< HEAD

int decrypt(unsigned char *e, const unsigned char *sk, const unsigned char *c);

#endif
=======
#define extract_preimage CRYPTO_NAMESPACE(extract_preimage)

int decrypt(unsigned char *e, const unsigned char *sk, const unsigned char *c);
int extract_preimage(unsigned char *preimage_out, const unsigned char *c, const unsigned char *sk, unsigned char *e);
#endif
>>>>>>> master
