
#include "utils.h"
#include "internal.h"

#include <stdarg.h>

static FILE *sx_log_handle_ = NULL;

void sx_set_log(FILE *io)
{
    sx_log_handle_ = io;
}

int sx_printf(const char *fmt, ...)
{
    va_list ap;
    int ret;

    va_start(ap, fmt);
    ret = vfprintf(stdout, fmt, ap);
    va_end(ap);

    if (sx_log_handle_ != NULL) {
        va_start(ap, fmt);
        vfprintf(sx_log_handle_, fmt, ap);
        va_end(ap);
    }

    fflush(stdout);
    fflush(sx_log_handle_);

    return ret;
}

/**
 * \fn void * sx_malloc(size_t size)
 *
 * \brief Allocate  memory
 *
 * \pars size    Number of memory blocks
 *
 * \return        Void pointer to the allocated memory
 *
 */
void * sx_malloc(size_t size)
{
    const size_t tsize = size;

    void *mem = NULL;

    if (tsize > 0) {
        mem = malloc(size);
    }

    if (mem == NULL) {
        sx_printf("### WARNING: Cannot allocate %.3"fFMTf" MB RAM!\n",
               (SX_FLT) tsize / 1048576);
    }

    return mem;
}

/**
 * \fn void * sx_calloc(size_t size, SX_INT type)
 *
 * \brief Allocate, initiate, and check memory
 *
 * \pars size    Number of memory blocks
 * \pars type    Size of memory blocks
 *
 * \return        Void pointer to the allocated memory
 *
 */
void * sx_calloc(size_t size, SX_INT type)
{
    const size_t tsize = size * type;

    void *mem = NULL;

    if (tsize > 0) {
        mem = calloc(size, type);
    }

    if (mem == NULL) {
        sx_printf("### WARNING: Cannot allocate %.3"fFMTf" MB RAM!\n",
               (SX_FLT) tsize / 1048576);
    }

    return mem;
}

/**
 * \fn void * sx_realloc (void * oldmem, size_t type)
 *
 * \brief Reallocate, initiate, and check memory
 *
 * \pars oldmem  Pointer to the existing mem block
 * \pars type    Size of memory blocks
 *
 * \return        Void pointer to the reallocated memory
 *
 */
void * sx_realloc(void *oldmem, size_t tsize)
{
    void *mem = NULL;

    if (tsize > 0) {
        mem = realloc(oldmem, tsize);
    }

    if (mem == NULL) {
        sx_printf("### WARNING: Cannot allocate %.3"fFMTf"MB RAM!\n",
               (SX_FLT) tsize / 1048576);
    }

    return mem;
}

/**
 * \fn void sx_free (void* mem)
 *
 * \brief Free up previous allocated memory body
 *
 * \pars mem   Pointer to the memory body need to be freed
 *
 * \return      NULL pointer
 *
 */
void sx_free(void *mem)
{
    if (mem) free(mem);
}

/**
 * \fn sx_get_time (SX_FLT *time)
 *
 * \brief Get system time
 *
 */
#if USE_UNIX

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <time.h>
#else
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#else
#include <time.h>
#endif
#endif

SX_FLT sx_get_time(void)
{
    struct timeval tv;
    double t;

    gettimeofday(&tv, (struct timezone *)0);
    t = tv.tv_sec + (double)tv.tv_usec * 1e-6;

    return t;
}

#else
#include <windows.h>

SX_FLT sx_get_time(void)
{
    LARGE_INTEGER timer;
    static LARGE_INTEGER fre;
    static int init = 0;
    double t;

    if (init != 1) {
        QueryPerformanceFrequency(&fre);
        init = 1;
    }

    QueryPerformanceCounter(&timer);

    t = timer.QuadPart * 1. / fre.QuadPart;

    return t;
}
#endif

#if USE_UNIX

#if HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif

SX_FLT sx_get_mem(void)
{
    struct rusage RU;
    SX_FLT mem_current;

    /* getrusage() */
    getrusage(RUSAGE_SELF, &RU);
    mem_current = RU.ru_maxrss / (double)1024.;

    return mem_current;
}

#else /* windows programming required */
SX_FLT sx_get_mem(void)
{
    return 0;
}

#endif

/**
 * \fn void sx_exit_on_errcode (const SX_INT status, const char *fctname)
 *
 * \brief Check error status and print out error messages before quit
 *
 * \pars status   Error status
 * \pars fctname  Function name where this routine is called
 *
 */
void sx_exit_on_errcode(const SX_INT status, const char *fctname)
{
    if (status >= 0) return;

    switch (status) {
        case ERROR_OPEN_FILE:
            sx_printf("### ERROR: %s -- Cannot open file!\n", fctname);
            break;

        case ERROR_WRONG_FILE:
            sx_printf("### ERROR: %s -- Wrong file format!\n", fctname);
            break;

        case ERROR_INPUT_PAR:
            sx_printf("### ERROR: %s -- Wrong input arguments!\n", fctname);
            break;

        case ERROR_ALLOC_MEM:
            sx_printf("### ERROR: %s -- Cannot allocate memory!\n", fctname);
            break;

        case ERROR_DATA_STRUCTURE:
            sx_printf("### ERROR: %s -- Data structure mismatch!\n", fctname);
            break;

        case ERROR_DATA_ZERODIAG:
            sx_printf("### ERROR: %s -- Matrix has zero diagonal entries!\n", fctname);
            break;

        case ERROR_DUMMY_VAR:
            sx_printf("### ERROR: %s -- Unexpected input argument!\n", fctname);
            break;

        case ERROR_AMG_INTERP_TYPE:
            sx_printf("### ERROR: %s -- Unknown AMG interpolation type!\n", fctname);
            break;

        case ERROR_AMG_COARSE_TYPE:
            sx_printf("### ERROR: %s -- Unknown AMG coarsening type!\n", fctname);
            break;

        case ERROR_AMG_SMOOTH_TYPE:
            sx_printf("### ERROR: %s -- Unknown AMG smoother type!\n", fctname);
            break;

        case ERROR_SOLVER_STAG:
            sx_printf("### ERROR: %s -- Solver stagnation error!\n", fctname);
            break;

        case ERROR_SOLVER_SOLSTAG:
            sx_printf("### ERROR: %s -- Solution is close to zero!\n", fctname);
            break;

        case ERROR_SOLVER_TOLSMALL:
            sx_printf("### ERROR: %s -- Tol is too small for the solver!\n", fctname);
            break;

        case ERROR_SOLVER_MAXIT:
            sx_printf("### ERROR: %s -- Max iteration number reached!\n", fctname);
            break;

        case ERROR_SOLVER_EXIT:
            sx_printf("### ERROR: %s -- Solver exited unexpected!\n", fctname);
            break;

        case ERROR_MISC:
            sx_printf("### ERROR: %s -- Unknown error occurred!\n", fctname);
            break;

        case ERROR_UNKNOWN:
            sx_printf("### ERROR: %s -- Function does not exit successfully!\n", fctname);
            break;

        default:
            break;
    }

    exit(status);
}

/**
 * \fn void sx_iarray_set (const SX_INT n, SX_INT *x, const SX_INT val)
 *
 * \brief Set initial value for an array to be x=val
 *
 * \pars n    Number of variables
 * \pars x    Pointer to the vector
 * \pars val  Initial value for the SX_FLT array
 *
 */
void sx_iarray_set(const SX_INT n, SX_INT *x, const SX_INT val)
{
    SX_INT i;

    assert(x != NULL);
    assert(n > 0);

    if (val == 0) {
        memset(x, 0, sizeof(SX_INT) * n);
    }
    else {
        for (i = 0; i < n; ++i) x[i] = val;
    }
}

/**
 * \fn void sx_iarray_cp (const SX_INT n, SX_INT *x, SX_INT *y) 
 *
 * \brief Copy an array to the other y=x
 *
 * \pars n    Number of variables
 * \pars x    Pointer to the original vector
 * \pars y    Pointer to the destination vector
 *
 */
void sx_iarray_cp(const SX_INT n, SX_INT *x, SX_INT *y)
{
    assert(x != NULL);
    assert(y != NULL);
    assert(n > 0);

    memcpy(y, x, n * sizeof(SX_INT));
}
