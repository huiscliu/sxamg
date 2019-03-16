
#ifndef SX_UTIL
#define SX_UTIL

#include "amg-type-def.h"

#ifdef __cplusplus
extern "C" {
#endif

void sx_set_log(FILE *io);
int sx_printf(const char *fmt, ...);

void * sx_mem_malloc(size_t size);
void * sx_mem_calloc(size_t size, SX_INT type);
void * sx_mem_realloc(void *oldmem, size_t tsize);
void sx_mem_free(void *mem);

SX_FLT sx_gettime(SX_FLT *time);
void sx_exit_on_errcode(const SX_INT status, const char *fctname);

#ifdef __cplusplus
}
#endif

#endif
