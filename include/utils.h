
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

#ifdef __cplusplus
}
#endif

#endif
