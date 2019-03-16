
#ifndef SX_INTERNAL_CONSTS
#define SX_INTERNAL_CONSTS

#include <assert.h>

/**
 * \brief Definition of logic type
 */
#define TRUE       1
#define FALSE      0

/**
 * \brief Definition of max, min, abs
 */
#define SX_MAX(a, b) (((a)>(b))?(a):(b))   /**< bigger one in a and b */
#define SX_MIN(a, b) (((a)<(b))?(a):(b))   /**< smaller one in a and b */
#define SX_ABS(a)    (((a)>=0.0)?(a):-(a)) /**< absolute value of a */

/**
 * \brief Type of vertices (DOFs) for coarsening
 */
#define UNPT                   -1  /**< Undetermined points */
#define FGPT                    0  /**< Fine grid points  */
#define CGPT                    1  /**< Coarse grid points */
#define ISPT                    2  /**< Isolated points */

#define BIGFLOAT             1e+20 /**< A large real number */
#define SMALLFLOAT           1e-20 /**< A small real number */
#define SMALLFLOAT2          1e-40 /**< An extremely small real number */
#define MAX_AMG_LVL            30  /**< Maximal AMG coarsening level */
#define MIN_CDOF               10  /**< Minimal number of coarsest variables */
#define MIN_CRATE             0.9  /**< Minimal coarsening ratio */
#define MAX_CRATE            20.0  /**< Maximal coarsening ratio */
#define MAX_RESTART            30  /**< Maximal restarting number */
#define MAX_STAG               20  /**< Maximal number of stagnation times */

/**
 * \brief Definition of iterative solver stopping criteria types
 */
typedef enum
{
    STOP_REL_RES        = 1,   /**< relative residual: ||r||/||r_0|| */
    STOP_REL_PRECRES    = 2,   /**< relative B-residual: ||r||_B/||b||_B */
    STOP_MOD_REL_RES    = 3,   /**< modified relative residual ||r||/||x|| */

} SX_STOP_TYPE;

typedef enum
{
    ERROR_OPEN_FILE          = -10,  /**< fail to open a file */
    ERROR_WRONG_FILE         = -11,  /**< input contains wrong format */
    ERROR_INPUT_PAR          = -12,  /**< wrong input argument */
    ERROR_MAT_SIZE           = -13,  /**< wrong problem size */
    ERROR_MISC               = -14,  /**< other error */

    ERROR_ALLOC_MEM          = -20,  /**< fail to allocate memory */
    ERROR_DATA_STRUCTURE     = -21,  /**< problem with data structures */
    ERROR_DATA_ZERODIAG      = -22,  /**< matrix has zero diagonal entries */
    ERROR_DUMMY_VAR          = -23,  /**< unexpected input data */

    ERROR_AMG_INTERP_TYPE    = -30,  /**< unknown interpolation type */
    ERROR_AMG_SMOOTH_TYPE    = -31,  /**< unknown smoother type */
    ERROR_AMG_COARSE_TYPE    = -32,  /**< unknown coarsening type */
    ERROR_AMG_COARSEING      = -33,  /**< coarsening step failed to complete */

    ERROR_SOLVER_STAG        = -42,  /**< solver stagnates */
    ERROR_SOLVER_SOLSTAG     = -43,  /**< solver's solution is too small */
    ERROR_SOLVER_TOLSMALL    = -44,  /**< solver's tolerance is too small */
    ERROR_SOLVER_MAXIT       = -48,  /**< maximal iteration number exceeded */
    ERROR_SOLVER_EXIT        = -49,  /**< solver does not quit successfully */

    ERROR_UNKNOWN            = -99,  /**< an unknown error type */

} SX_ERROR_CODE;

#ifdef __cplusplus
extern "C" {
#endif

void sx_print_itinfo(const SX_INT ptrlvl, const SX_INT stop_type, const SX_INT iter,
        const SX_FLT relres, const SX_FLT absres, const SX_FLT factor);

void sx_iarray_cp(const SX_INT n, SX_INT *x, SX_INT *y);
void sx_iarray_set(const SX_INT n, SX_INT *x, const SX_INT val);

/* matrix */
SX_IMAT sx_imat_trans(SX_IMAT *A);
void sx_imat_cp(SX_IMAT *A, SX_IMAT *B);
void sx_imat_destroy(SX_IMAT *A);
SX_IMAT sx_imat_struct_create(SX_INT nrow, SX_INT ncol, SX_INT nnz);

/* vector */
SX_IVEC sx_ivec_create(SX_INT m);
void sx_ivec_destroy(SX_IVEC *u);
void sx_ivec_set(SX_IVEC *u, SX_INT val);

#ifdef __cplusplus
}
#endif

#endif
