
#ifndef SX_AMG_SETUP
#define SX_AMG_SETUP

#include "amg-type-def.h"
#include "vec.h"
#include "op-coarsen.h"
#include "op-inter.h"
#include "amg-utils.h"

#ifdef __cplusplus
extern "C" {
#endif

void sx_amg_setup(SX_AMG *mg, SX_MAT *A, SX_AMG_PARS *pars);

#ifdef __cplusplus
}
#endif

#endif
