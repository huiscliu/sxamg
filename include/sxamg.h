
#ifndef SX_SXAMG
#define SX_SXAMG

#include "amg-type-def.h"
#include "utils.h"
#include "op-blas.h"
#include "amg-utils.h"
#include "amg-setup.h"
#include "amg-cycle.h"

#ifdef __cplusplus
extern "C" {
#endif

SX_RTN sx_solver_amg_solve(SX_AMG *mg, SX_VEC *x, SX_VEC *b);
SX_RTN sx_solver_amg(SX_MAT *A, SX_VEC *x, SX_VEC *b, SX_AMG_PARS *pars);

#ifdef __cplusplus
}
#endif

#endif
