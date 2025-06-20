#ifndef PARAMS_H
#define PARAMS_H

#include "namespace.h"

/*
#define GFBITS 13
#define SYS_N 8192
#define SYS_T 128
*/


#define GFBITS 12
#define SYS_N 3488
#define SYS_T 64

#define COND_BYTES ((1 << (GFBITS-4))*(2*GFBITS - 1))
#define IRR_BYTES (SYS_T * 2)

#define PK_NROWS (SYS_T*GFBITS) // equivalent to k
#define PK_NCOLS (SYS_N - PK_NROWS)

#define PK_ROW_BYTES ((PK_NCOLS + 7)/8)

#define SYND_BYTES ((PK_NROWS + 7)/8)

#define GFMASK ((1 << GFBITS) - 1)

#endif


// PK_BYTES = SYS_T*GFBITS*(SYS_N-SYS_T*GFBITS)/8