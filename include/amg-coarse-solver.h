
#ifndef SX_CS_SOLVER
#define SX_CS_SOLVER

#include "amg-type-def.h"
#include "utils.h"
#include "op-blas.h"

#ifdef __cplusplus
extern "C" {
#endif

void sx_amg_coarest_solve(SX_AMG_COMP *cg, const SX_FLT ctol, SX_CSSOLVE_TYPE cssolve_type, const SX_INT verb);

#ifdef __cplusplus
}
#endif

#endif
