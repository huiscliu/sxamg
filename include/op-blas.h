
#ifndef SX_OP_BLAS
#define SX_OP_BLAS

#include "amg-type-def.h"
#include "utils.h"
#include "mat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* array */
SX_FLOAT sx_blas_array_norminf(SX_INT n, const SX_FLOAT *x);
SX_FLOAT sx_blas_array_norm2(SX_INT n, const SX_FLOAT *x);
SX_FLOAT sx_blas_array_norm1(SX_INT n, const SX_FLOAT *x);
SX_FLOAT sx_blas_array_dot(SX_INT n, const SX_FLOAT *x, const SX_FLOAT *y);

void sx_blas_array_axpby(SX_INT n, SX_FLOAT a, const SX_FLOAT *x, SX_FLOAT b, SX_FLOAT *y);
void sx_blas_array_axpyz(SX_INT n, SX_FLOAT a, const SX_FLOAT *x, const SX_FLOAT *y, SX_FLOAT *z);
void sx_blas_array_axpy(SX_INT n, SX_FLOAT a, const SX_FLOAT *x, SX_FLOAT *y);
void sx_blas_array_ax(SX_INT n, SX_FLOAT a, SX_FLOAT *x);
void sx_blas_array_cp(SX_INT n, const SX_FLOAT *x, SX_FLOAT *y);
void sx_blas_array_set(SX_INT n, SX_FLOAT *x, SX_FLOAT val);

/* vector */
SX_FLOAT sx_blas_vec_norm2(const SX_VEC *x);
SX_FLOAT sx_blas_vec_dot(const SX_VEC *x, const SX_VEC *y);

void sx_blas_vec_axpyz(SX_FLOAT a, const SX_VEC *x, const SX_VEC *y, SX_VEC *z);
void sx_blas_vec_axpbyz(SX_FLOAT a, const SX_VEC *x, SX_FLOAT b, const SX_VEC *y, SX_VEC *z);

void sx_blas_vec_axpy(SX_FLOAT a, const SX_VEC *x, SX_VEC *y);
void sx_blas_vec_axpby(SX_FLOAT a, const SX_VEC *x, SX_FLOAT b, SX_VEC *y);

/* mat-vec */
void sx_blas_mat_amxpy(SX_FLOAT a, const SX_MAT *A, const SX_VEC *x, SX_VEC *y);
void sx_blas_mat_mxy(const SX_MAT *A, const SX_VEC *x, SX_VEC *y);

/* mat-mat */
SX_MAT sx_blas_mat_rap(const SX_MAT *R, const SX_MAT *A, const SX_MAT *P);

#ifdef __cplusplus
}
#endif

#endif
