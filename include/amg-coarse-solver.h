
#ifndef SX_CS_SOLVER
#define SX_CS_SOLVER

#include "amg-type-def.h"
#include "utils.h"
#include "op-blas.h"

#ifdef __cplusplus
extern "C" {
#endif

void sx_amg_coarest_solve(SX_MAT *A, SX_VEC *b, SX_VEC *x, const SX_FLT ctol, const SX_INT verb);

#ifdef __cplusplus
}
#endif

#endif
