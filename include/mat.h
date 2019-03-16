
#ifndef SX_MAT_H
#define SX_MAT_H

#include "amg-type-def.h"
#include "vec.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

SX_MAT sx_mat_struct_create(SX_INT nrow, SX_INT ncol, SX_INT nnz);
SX_MAT sx_mat_create(SX_INT nrow, SX_INT ncol, SX_INT *Ap, SX_INT *Aj, SX_FLT *Ax);

void sx_mat_destroy(SX_MAT *A);

SX_MAT sx_mat_trans(SX_MAT *A);
void sx_mat_sort(SX_MAT *A);

void sx_mat_cp(SX_MAT *A, SX_MAT *B);

SX_VEC sx_mat_get_diag(SX_MAT *A, SX_INT n);

#ifdef __cplusplus
}
#endif

#endif
