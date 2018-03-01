
#ifndef SX_OP_COARSEN
#define SX_OP_COARSEN

#include "amg-type-def.h"
#include "utils.h"
#include "mat.h"
#include "op-blas.h"
#include "vec.h"

#ifdef __cplusplus
extern "C" {
#endif

SX_INT sx_amg_coarsen(SX_MAT *A, SX_IVEC *vertices, SX_MAT *P, SX_IMAT *S, SX_AMG_PARS *pars);

#ifdef __cplusplus
}
#endif

#endif
