
#ifndef SX_OP_BLAS
#define SX_OP_BLAS

#include "amg-type-def.h"
#include "utils.h"
#include "mat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* array */
SX_FLT sx_blas_array_norminf(SX_INT n, const SX_FLT *x);
SX_FLT sx_blas_array_norm2(SX_INT n, const SX_FLT *x);
SX_FLT sx_blas_array_norm1(SX_INT n, const SX_FLT *x);
SX_FLT sx_blas_array_dot(SX_INT n, const SX_FLT *x, const SX_FLT *y);

void sx_blas_array_axpby(SX_INT n, SX_FLT a, const SX_FLT *x, SX_FLT b, SX_FLT *y);
void sx_blas_array_axpyz(SX_INT n, SX_FLT a, const SX_FLT *x, const SX_FLT *y, SX_FLT *z);
void sx_blas_array_axpy(SX_INT n, SX_FLT a, const SX_FLT *x, SX_FLT *y);
void sx_blas_array_ax(SX_INT n, SX_FLT a, SX_FLT *x);
void sx_blas_array_cp(SX_INT n, const SX_FLT *x, SX_FLT *y);
void sx_blas_array_set(SX_INT n, SX_FLT *x, SX_FLT val);

/* vector */
SX_FLT sx_blas_vec_norm2(const SX_VEC *x);
SX_FLT sx_blas_vec_dot(const SX_VEC *x, const SX_VEC *y);

void sx_blas_vec_axpyz(SX_FLT a, const SX_VEC *x, const SX_VEC *y, SX_VEC *z);
void sx_blas_vec_axpbyz(SX_FLT a, const SX_VEC *x, SX_FLT b, const SX_VEC *y, SX_VEC *z);

void sx_blas_vec_axpy(SX_FLT a, const SX_VEC *x, SX_VEC *y);
void sx_blas_vec_axpby(SX_FLT a, const SX_VEC *x, SX_FLT b, SX_VEC *y);

void sx_blas_vec_copy(const SX_VEC *x, SX_VEC *y);
void sx_blas_vec_set(SX_VEC *x, SX_FLT val);

/* mat-vec */
void sx_blas_mv_amxpy(SX_FLT a, const SX_MAT *A, const SX_VEC *x, SX_VEC *y);
void sx_blas_mv_amxpby(SX_FLT a, const SX_MAT *A, const SX_VEC *x, SX_FLT b, SX_VEC *y);
void sx_blas_mv_amxpbyz(SX_FLT a, const SX_MAT *A, const SX_VEC *x, SX_FLT b, SX_VEC *y,  SX_VEC *z);

void sx_blas_mv_mxy(const SX_MAT *A, const SX_VEC *x, SX_VEC *y);

/* mat-mat: C = RAP */
SX_MAT sx_blas_mat_rap(const SX_MAT *R, const SX_MAT *A, const SX_MAT *P);

/* mat-mat: C = AB */
void sx_blas_mat_mxm (const SX_MAT  *A, const SX_MAT  *B, SX_MAT *C);

#ifdef __cplusplus
}
#endif

#endif
