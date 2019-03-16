
#ifndef SX_SOLVER_VEC
#define SX_SOLVER_VEC

#include "amg-type-def.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

SX_VEC sx_vec_create(SX_INT m);
void sx_vec_destroy(SX_VEC *u);

void sx_vec_set_value(SX_VEC *x, SX_FLT val);
void sx_vec_cp(const SX_VEC *src, SX_VEC *des);
SX_INT sx_vec_get_size(const SX_VEC *v);

void sx_vec_set_entry(SX_VEC *x, SX_INT index, SX_FLT val);
SX_FLT sx_vec_get_entry(const SX_VEC *x, SX_INT index);

#ifdef __cplusplus
}
#endif

#endif
