
#ifndef SX_UTIL
#define SX_UTIL

#include "amg-type-def.h"

#ifdef __cplusplus
extern "C" {
#endif

void sx_set_log(FILE *io);
int sx_printf(const char *fmt, ...);

void * sx_malloc(size_t size);
void * sx_calloc(size_t size, SX_INT type);
void * sx_realloc(void *oldmem, size_t tsize);
void sx_free(void *mem);

SX_FLT sx_get_time(void);
SX_FLT sx_get_mem(void);

void sx_exit_on_errcode(const SX_INT status, const char *fctname);

/* LU factorization */
int sx_mat_dense_lu(int n, SX_FLT *A, int pvt[]);

/* solves L*U*X = B  */
void sx_mat_dense_sv(int n, SX_FLT *A, int pvt[], int m, SX_FLT *B);

/* solves AX = B, returns 1 if successful and 0 if A is singular */
int sx_mat_dense_solve(int n, int m, SX_FLT *A, SX_FLT *B);

#ifdef __cplusplus
}
#endif

#endif
