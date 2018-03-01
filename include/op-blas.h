
#ifndef SX_OP_BLAS
#define SX_OP_BLAS

#include "amg-type-def.h"
#include "utils.h"
#include "mat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* array */
SX_FLOAT sx_blas_array_norminf(const SX_INT n, const SX_FLOAT *x);
SX_FLOAT sx_blas_array_norm2(const SX_INT n, const SX_FLOAT *x);
SX_FLOAT sx_blas_array_norm1(const SX_INT n, const SX_FLOAT *x);
SX_FLOAT sx_blas_array_dot(const SX_INT n, const SX_FLOAT *x, const SX_FLOAT *y);

void sx_blas_array_axpby(const SX_INT n, const SX_FLOAT a, SX_FLOAT *x, const SX_FLOAT b, SX_FLOAT *y);
void sx_blas_array_axpyz(const SX_INT n, const SX_FLOAT a, SX_FLOAT *x, SX_FLOAT *y, SX_FLOAT *z);
void sx_blas_array_axpy(const SX_INT n, const SX_FLOAT a, SX_FLOAT *x, SX_FLOAT *y);
void sx_blas_array_ax(const SX_INT n, const SX_FLOAT a, SX_FLOAT *x);
void sx_blas_array_cp(const SX_INT n, SX_FLOAT *x, SX_FLOAT *y);
void sx_blas_array_set(const SX_INT n, SX_FLOAT *x, const SX_FLOAT val);

/* vector */
SX_FLOAT sx_blas_vec_norm2(SX_VEC *x);
SX_FLOAT sx_blas_vec_dot(SX_VEC *x, SX_VEC *y);

void sx_blas_vec_axpyz(const SX_FLOAT a, SX_VEC *x, SX_VEC *y, SX_VEC *z);
void sx_blas_vec_axpy(const SX_FLOAT a, SX_VEC *x, SX_VEC *y);

/* mat-vec */
void sx_blas_mat_amxpy(const SX_FLOAT alpha, SX_MAT *A, SX_VEC *x, SX_VEC *y);
void sx_blas_mat_mxy(SX_MAT *A, SX_VEC *x, SX_VEC *y);

/* mat-mat */
SX_MAT sx_blas_mat_rap(SX_MAT *R, SX_MAT *A, SX_MAT *P);

#ifdef __cplusplus
}
#endif

#endif
