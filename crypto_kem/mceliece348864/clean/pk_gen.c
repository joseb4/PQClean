/*
  This file is for public-key generation
*/

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "controlbits.h"
#include "benes.h"
#include "crypto_declassify.h"
#include "crypto_uint64.h"
#include "params.h"
#include "pk_gen.h"
#include "root.h"
#include "uint64_sort.h"
#include "util.h"

static crypto_uint64 uint64_is_equal_declassify(uint64_t t, uint64_t u) {
    crypto_uint64 mask = crypto_uint64_equal_mask(t, u);
    crypto_declassify(&mask, sizeof mask);
    return mask;
}

static crypto_uint64 uint64_is_zero_declassify(uint64_t t) {
    crypto_uint64 mask = crypto_uint64_zero_mask(t);
    crypto_declassify(&mask, sizeof mask);
    return mask;
}

/* input: secret key sk */
/* output: public key pk */
int pk_gen(unsigned char *pk, unsigned char *sk, const uint32_t *perm, int16_t *pi) {
    int i, j, k;
    int row, c;

    uint64_t buf[ 1 << GFBITS ];

    unsigned char mat[ PK_NROWS ][ SYS_N / 8 ];
    unsigned char mask;
    unsigned char b;

    gf g[ SYS_T + 1 ]; // Goppa polynomial
    gf L[ SYS_N ]; // support
    gf inv[ SYS_N ];

    //

    g[ SYS_T ] = 1;

    for (i = 0; i < SYS_T; i++) {
        g[i] = load_gf(sk);
        sk += 2;
    }

    for (i = 0; i < (1 << GFBITS); i++) {
        buf[i] = perm[i];
        buf[i] <<= 31;
        buf[i] |= i;
    }

    uint64_sort(buf, 1 << GFBITS);

    for (i = 1; i < (1 << GFBITS); i++) {
        if (uint64_is_equal_declassify(buf[i - 1] >> 31, buf[i] >> 31)) {
            return -1;
        }
    }

    for (i = 0; i < (1 << GFBITS); i++) {
        pi[i] = buf[i] & GFMASK;
    }
    for (i = 0; i < SYS_N;         i++) {
        L[i] = bitrev(pi[i]);
    }

    // filling the matrix

    root(inv, g, L);

    for (i = 0; i < SYS_N; i++) {
        inv[i] = gf_inv(inv[i]);
    }

    for (i = 0; i < PK_NROWS; i++) {
        for (j = 0; j < SYS_N / 8; j++) {
            mat[i][j] = 0;
        }
    }

    for (i = 0; i < SYS_T; i++) {
        for (j = 0; j < SYS_N; j += 8) {
            for (k = 0; k < GFBITS;  k++) {
                b  = (inv[j + 7] >> k) & 1;
                b <<= 1;
                b |= (inv[j + 6] >> k) & 1;
                b <<= 1;
                b |= (inv[j + 5] >> k) & 1;
                b <<= 1;
                b |= (inv[j + 4] >> k) & 1;
                b <<= 1;
                b |= (inv[j + 3] >> k) & 1;
                b <<= 1;
                b |= (inv[j + 2] >> k) & 1;
                b <<= 1;
                b |= (inv[j + 1] >> k) & 1;
                b <<= 1;
                b |= (inv[j + 0] >> k) & 1;

                mat[ i * GFBITS + k ][ j / 8 ] = b;
            }
        }

        for (j = 0; j < SYS_N; j++) {
            inv[j] = gf_mul(inv[j], L[j]);
        }

    }

    // gaussian elimination

    for (i = 0; i < (PK_NROWS + 7) / 8; i++) {
        for (j = 0; j < 8; j++) {
            row = i * 8 + j;

            if (row >= PK_NROWS) {
                break;
            }

            for (k = row + 1; k < PK_NROWS; k++) {
                mask = mat[ row ][ i ] ^ mat[ k ][ i ];
                mask >>= j;
                mask &= 1;
                mask = -mask;

                for (c = 0; c < SYS_N / 8; c++) {
                    mat[ row ][ c ] ^= mat[ k ][ c ] & mask;
                }
            }

            if ( uint64_is_zero_declassify((mat[ row ][ i ] >> j) & 1) ) { // return if not systematic
                //printf("Error: row %d is not systematic\n", row);
                return -1;
            }

            for (k = 0; k < PK_NROWS; k++) {
                if (k != row) {
                    mask = mat[ k ][ i ] >> j;
                    mask &= 1;
                    mask = -mask;

                    for (c = 0; c < SYS_N / 8; c++) {
                        mat[ k ][ c ] ^= mat[ row ][ c ] & mask;
                    }
                }
            }
        }
    }

    for (i = 0; i < PK_NROWS; i++) {
        memcpy(pk + i * PK_ROW_BYTES, mat[i] + PK_NROWS / 8, PK_ROW_BYTES);
    }

    return 0;
}



/*
 * Construye la matriz generadora G en forma sistemática a partir de la clave pública.
 *
 * La clave pública pk contiene la matriz A (de dimensión 
 *   PK_NROWS x ((SYS_N - PK_NROWS) bits)
 * empaquetada en bytes (cada fila ocupa PK_ROW_BYTES bytes).
 *
 * La matriz generadora G se construye como:
 *      G = [ -A^T | I ]
 * de dimensión: (SYS_N - PK_NROWS) x SYS_N bits,
 * es decir, (SYS_N - PK_NROWS) filas, y cada fila se almacena en SYS_N/8 bytes.
 *
 * Parámetros:
 *   G  : salida, matriz generadora, definida como: 
 *        unsigned char G[SYS_N - PK_NROWS][SYS_N/8]
 *   pk : entrada, clave pública que contiene la matriz A, de dimensión:
 *        unsigned char pk[PK_NROWS * PK_ROW_BYTES]
 */
void build_G_matrix(unsigned char G[][SYS_N / 8], const unsigned char *pk) {
    int i, j;
    // Número de filas de G (la dimensión k del código) es:
    int num_rows_G = SYS_N - PK_NROWS;  // k = n - (n - k) = SYS_N - PK_NROWS
    int total_bytes = SYS_N / 8;         // número de bytes por fila de G
    
    // Variable temporal para construir cada fila de G
    unsigned char row[SYS_N / 8];
    
    /*
     * Para cada fila i de G (0 <= i < num_rows_G), se construye:
     *  - La parte izquierda: A^T, es decir, para cada columna j (0 <= j < PK_NROWS),
     *    se extrae el bit de la fila j de A, en la posición i.
     *  - La parte derecha: la matriz identidad, se coloca un 1 en la posición (PK_NROWS + i)
     *    de la fila.
     */
    for (i = 0; i < num_rows_G; i++) {
        // Inicializar la fila a cero.
        for (j = 0; j < total_bytes; j++) {
            row[j] = 0;
        }
        
        // Construir la parte izquierda: para cada j, la entrada (j, i) de A^T es el bit (i) de la fila j de A.
        // La matriz A está almacenada en pk de forma consecutiva: 
        //   A[j] se encuentra en pk + j * PK_ROW_BYTES, y contiene PK_NCOLS bits, donde PK_NCOLS = SYS_N - PK_NROWS.
        for (j = 0; j < PK_NROWS; j++) {
            const unsigned char *A_row = pk + j * PK_ROW_BYTES;
            // Extraer el bit i de A_row.
            int bit = (A_row[i / 8] >> (i % 8)) & 1;
            if (bit)
                row[j / 8] |= (1 << (j % 8));
        }
        
        // Construir la parte derecha: identidad de dimensión (SYS_N - PK_NROWS) bits.
        // En la fila i, se coloca un 1 en la posición (PK_NROWS + i).
        int pos = PK_NROWS + i;
        row[pos / 8] |= (1 << (pos % 8));
        
        // Copiar la fila construida en G[i].
        memcpy(G[i], row, total_bytes);
    }
}
