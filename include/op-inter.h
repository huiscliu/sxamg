
#ifndef SX_OP_INTER
#define SX_OP_INTER

#include "amg-type-def.h"
#include "utils.h"
#include "op-blas.h"

#ifdef __cplusplus
extern "C" {
#endif

void sx_amg_interp_trunc(SX_MAT *P, SX_AMG_PARS *pars);

void sx_amg_interp(SX_MAT *A, SX_IVEC *vertices, SX_MAT *P, SX_IMAT *S, SX_AMG_PARS *pars);

#ifdef __cplusplus
}
#endif

#endif
