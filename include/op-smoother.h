
#ifndef SX_OP_SMOOTH
#define SX_OP_SMOOTH

#include "amg-type-def.h"
#include "utils.h"
#include "op-blas.h"

#ifdef __cplusplus
extern "C" {
#endif

void sx_amg_smoother_post(SX_SMTR *s);
void sx_amg_smoother_pre(SX_SMTR *s);

#ifdef __cplusplus
}
#endif

#endif
