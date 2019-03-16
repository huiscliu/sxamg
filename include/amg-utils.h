
#ifndef SX_AMG_UTIL
#define SX_AMG_UTIL

#include "amg-type-def.h"
#include "utils.h"
#include "mat.h"
#include "amg-coarse-solver.h"
#include "op-smoother.h"

#ifdef __cplusplus
extern "C" {
#endif

SX_AMG sx_amg_data_create(SX_AMG_PARS *pars);
void sx_amg_data_destroy(SX_AMG *mg);

void sx_amg_pars_print(SX_AMG_PARS *pars);
void sx_amg_pars_init(SX_AMG_PARS *pars);
void sx_amg_complexity_print(SX_AMG *mg);

#ifdef __cplusplus
}
#endif

#endif
