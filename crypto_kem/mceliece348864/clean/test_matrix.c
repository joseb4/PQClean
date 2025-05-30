#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include "api.h"
#include "operations.h"
#include "params.h"

#include <stdio.h>




// Función para manejar errores
void handle_error(const char *msg) {
    fprintf(stderr, "Error: %s\n", msg);
    exit(EXIT_FAILURE);
}
/* input: clave pública pk */
/* output: Matriz H con ( I || A ) en bytes */
void build_H_matrix(unsigned char H[][SYS_N / 8], const unsigned char *pk) {
    unsigned char row[SYS_N / 8];
    const unsigned char *pk_ptr = pk;

    int i, j;

    for (i = 0; i < PK_NROWS; i++) {
        // Inicializar la fila en 0
        for (j = 0; j < SYS_N / 8; j++) {
            row[j] = 0;
        }

        // Copiar la parte A (desde la clave pública)
        for (j = 0; j < PK_ROW_BYTES; j++) {
            row[ SYS_N / 8 - PK_ROW_BYTES + j ] = pk_ptr[j];
        }

        // Insertar el bit correspondiente a la matriz identidad en la parte izquierda
        // row[i / 8] |= 1 << (i % 8);
        row[i / 8] |= 0x80 >> (i % 8);


        // Copiar esta fila completa en la matriz H de salida
        memcpy(H[i], row, SYS_N / 8);

        // Avanzar el puntero en la clave pública para la siguiente fila de A
        pk_ptr += PK_ROW_BYTES;
    }
}

/* input: clave pública pk */
/* output: Matriz G_t con ( A; I ) en bytes, es decir G_t = (A^t || I)^t = [A; I] */
void build_Gt_matrix(unsigned char Gt[][ (SYS_N - PK_NROWS) / 8 ], const unsigned char *pk) {
    int i, j;
    const unsigned char *pk_ptr = pk;

    /* 
     * Parte 1: Copiar la matriz A
     * La matriz A está representada en pk y ocupa PK_NROWS renglones, 
     * cada uno con PK_ROW_BYTES bytes.
     */
    for (i = 0; i < PK_NROWS; i++) {
        memcpy(Gt[i], pk_ptr, PK_ROW_BYTES);
        pk_ptr += PK_ROW_BYTES;
    }

    /*
     * Parte 2: Concatenar la matriz identidad por abajo.
     * Se generan SYS_N - PK_NROWS renglones (la dimensión de I en la forma sistemática de G_t),
     * donde cada renglón tiene (SYS_N - PK_NROWS) bits (o (SYS_N - PK_NROWS)/8 bytes).
     */
    int identity_rows = SYS_N - PK_NROWS;
    for (i = 0; i < identity_rows; i++) {
        /* Inicializar el renglón a 0 */
        for (j = 0; j < (SYS_N - PK_NROWS) / 8; j++) {
            Gt[PK_NROWS + i][j] = 0;
        }
        /* Activar el bit correspondiente a la identidad:
         * En la fila i de la parte identidad, se coloca un 1 en la posición i.
         * Esto se logra accediendo al byte i/8 y activando el bit (i % 8).
         */
        Gt[PK_NROWS + i][ i / 8 ] |= (0x80 >> (i % 8));

    }
}

void matrix_mult(int m, int n, int p, 
                 unsigned char A[m][n], 
                 unsigned char B[n][p], 
                 unsigned char result[m][p]) {
    int i, j, k;
    // Se recorre cada fila de A y cada columna de B
    for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
            result[i][j] = 0; // Inicializar la celda del resultado
            for (k = 0; k < n; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}



// Función para transponer una matriz.
// Parámetros:
//   n_rows: número de filas de la matriz original.
//   n_cols: número de columnas de la matriz original.
//   matriz: matriz original de tamaño n_rows x n_cols.
//   matrizT: matriz destino donde se almacenará la transpuesta (tamaño n_cols x n_rows).
void traspose(int n_rows, int n_cols, 
    unsigned char matriz[n_rows][n_cols], 
    unsigned char matrizT[n_cols][n_rows]) {
    int i, j;
    for (i = 0; i < n_rows; i++) {
        for (j = 0; j < n_cols; j++) {
            matrizT[j][i] = matriz[i][j];
        }
    }
}


/*
 * Construye la matriz generadora G en forma sistemática a partir de la clave pública.
 *
 * La clave pública pk contiene la matriz A (de dimensión 
 *   PK_NROWS x ((SYS_N - PK_NROWS) bits)
 * empaquetada en bytes (cada fila ocupa PK_ROW_BYTES bytes).
 *
 * La matriz generadora G se construye como:
 *      G = [ A^T | I ]
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
void build_G_transpose(unsigned char G_t[][ (SYS_N - PK_NROWS + 7) / 8 ], const unsigned char *pk) {
    int num_rows_G = SYS_N - PK_NROWS;  // Número de filas de G.
    int out_cols = (num_rows_G + 7) / 8;  // Número de bytes por fila en G^T.
    int r, c;

    // Recorrer cada fila "r" de G^T, que corresponde a la columna "r" de G.
    for (r = 0; r < SYS_N; r++) {
        // Inicializamos la fila r de G_t a cero.
        for (int j = 0; j < out_cols; j++) {
            G_t[r][j] = 0;
        }
        // Ahora, para cada columna "c" (0 <= c < num_rows_G) de G^T,
        // se asigna: G_t[r][c] = G[c][r].
        for (c = 0; c < num_rows_G; c++) {
            int bit;
            if (r < PK_NROWS) {
                // Caso 1: r corresponde a la parte de A^T en G.
                // G[c][r] es el bit "c" de la fila "r" de A.
                const unsigned char *pk_row = pk + r * PK_ROW_BYTES;
                bit = (pk_row[c / 8] >> (c % 8)) & 1;
            } else {
                // Caso 2: r corresponde a la parte identidad de G.
                // Se tiene que r = PK_NROWS + t, donde t = r - PK_NROWS.
                // En G, la identidad se coloca en la posición PK_NROWS + c para la fila c,
                // por lo que G[c][r] = 1 si y solo si r == PK_NROWS + c.
                int t = r - PK_NROWS;
                bit = (t == c) ? 1 : 0;
            }
            if (bit)
                G_t[r][c / 8] |= (1 << (c % 8));
        }
    }
}


void print_hex_matrix(const char *name, const unsigned char *matrix, int rows, int cols) {
    printf("Matrix %s (%d x %d bytes):\n", name, rows, cols);
    for (int i = 0; i < rows; i++) {
        printf("[%3d] ", i);
        for (int j = 0; j < cols; j++) {
            printf("%02X ", matrix[i * cols + j]);
        }
        printf("\n");
    }
}

/*
 * Multiplica dos matrices binarias sobre GF(2).
 *
 * Parámetros:
 *   m: número de filas de la matriz A (y de C).
 *   n: número de columnas (en bits) de A y número de filas (en bits) de B.
 *   p: número de columnas (en bits) de la matriz B (y de C).
 *   A: matriz A, de dimensión m x n bits, empaquetada en filas de ((n+7)/8) bytes.
 *   B: matriz B, de dimensión n x p bits, empaquetada en filas de ((p+7)/8) bytes.
 *   C: matriz resultado (producto A x B), de dimensión m x p bits, empaquetada en filas de ((p+7)/8) bytes.
 *
 * Nota: Se asume C99 o posterior, para permitir VLA (Variable Length Arrays)
 * en los parámetros de la función.
 */
void multiply_matrices_GF2(int m, int n, int p,
                           const unsigned char A[][ (n + 7) / 8 ],
                           const unsigned char B[][ (p + 7) / 8 ],
                           unsigned char C[][ (p + 7) / 8 ]) {
    int i, j, k;
    int p_bytes = (p + 7) / 8;
    
    // Inicializar la matriz C a cero.
    for (i = 0; i < m; i++) {
        for (j = 0; j < p_bytes; j++) {
            C[i][j] = 0;
        }
    }
    
    // Multiplicación de las matrices en GF(2):
    // Para cada entrada C[i][j] se calcula la suma (XOR) de A[i,k] * B[k,j] para k=0..n-1.
    for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
            int dot = 0;
            for (k = 0; k < n; k++) {
                // Extraemos el k-ésimo bit de la fila i de A.
                int bitA = (A[i][k / 8] >> (k % 8)) & 1;
                // Extraemos el k-ésimo bit (de la fila k) de B en la columna j.
                int bitB = (B[k][j / 8] >> (j % 8)) & 1;
                // La operación AND y la suma modulo 2 (XOR).
                dot ^= (bitA & bitB);
            }
            // Si el resultado del producto interno (dot) es 1, se activa el bit correspondiente
            if (dot) {
                C[i][j / 8] |= (1 << (j % 8));
            }
        }
    }
}

#include <stdio.h>
#include <stdlib.h>

void save_matrix_raw(const char *filename, const unsigned char *matrix, size_t rows, size_t cols) {
    FILE *f = fopen(filename, "wb");
    if (!f) {
        perror("Error abriendo el archivo para escritura");
        exit(EXIT_FAILURE);
    }

    size_t written = fwrite(matrix, 1, rows * cols, f);
    if (written != rows * cols) {
        fprintf(stderr, "Error al escribir la matriz completa en el archivo\n");
        exit(EXIT_FAILURE);
    }

    fclose(f);
}
// Función que multiplica la matriz binaria 'matriz_a' por el vector 'vector',
// todo en binario (mod 2), y guarda el resultado en 'resultado'.
void multiplicar_mod2(unsigned char matriz_a[SYS_N][(SYS_N - PK_NROWS) / 8],
                      unsigned char vector_col[(SYS_N - PK_NROWS)][1],
                      unsigned char resultado[SYS_N][1]) 
{
    for (int i = 0; i < SYS_N; i++) {
        unsigned char acum = 0;
        for (int j = 0; j < (SYS_N - PK_NROWS); j++) {
            // En binario, la multiplicación es (matriz_a[i][j] & vector_col[j][0])
            // y la suma mod 2 es XOR (acum ^= ...)
            acum ^= (matriz_a[i][j] & vector_col[j][0]);
        }
        resultado[i][0] = acum;  // resultado binario para la fila i
    }
}
// Imprime un vector columna binario (unsigned char) en formato hexadecimal
void print_col_hex(unsigned char vector_col[][1], int nfilas) {
    for (int i = 0; i < nfilas; i++) {
        // %02X -> imprime en hexadecimal (base 16), con dos dígitos y en mayúsculas
        // Ej: si vector_col[i][0] es 5, se imprime "05"; si es 15, "0F".
        printf("%02X", vector_col[i][0]);
    }
}
int main(){
    srand((unsigned) time(NULL));
    unsigned char pk[PQCLEAN_MCELIECE6688128F_CLEAN_CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[PQCLEAN_MCELIECE6688128F_CLEAN_CRYPTO_SECRETKEYBYTES];

    if (crypto_kem_keypair(pk, sk) != 0) {
        printf("Error generando claves de McEliece");
    }

    int k = SYS_N - PK_NROWS;
    unsigned char mensaje[k][1];
    unsigned char byte = 0xAB; // 10101011
    

    for (int i = 0; i < k; i++) {
        // Extrae el bit i del byte, desde el bit más significativo (7) al menos (0)
        mensaje[i][0] = (byte >> (7 - i)) & 1;
    }
    print_col_hex(mensaje, k);


    unsigned char H_matrix[PK_NROWS][SYS_N / 8];
    build_H_matrix(H_matrix, pk);

    
    unsigned char G_transpose[SYS_N] [(SYS_N - PK_NROWS) / 8];
    build_Gt_matrix(G_transpose, pk);

    // save_matrix_raw("H_matrix.raw", (const unsigned char *)H_matrix, PK_NROWS, SYS_N/8);
    // save_matrix_raw("G_transpose.raw", (const unsigned char *)G_transpose, SYS_N, (SYS_N - PK_NROWS)/8);

    unsigned char resultado[SYS_N][1];
    multiplicar_mod2(G_transpose, mensaje, resultado);
    // print_col_hex(resultado, SYS_N);
    // print_hex_matrix("H", (const unsigned char *)H_matrix, 10/* PK_NROWS */,SYS_N/8);
    // print_hex_matrix("G_transpuesta", (const unsigned char *)G_transpose, SYS_N, (SYS_N-PK_NROWS)/8);
    // unsigned char G_matrix[SYS_N - PK_NROWS][SYS_N / 8];
    // build_G_matrix(G_matrix, pk);
    // Recuerda liberar la memoria con free() cuando ya no la necesites.


    
    return 0;
}