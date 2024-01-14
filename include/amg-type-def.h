
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifndef SX_TYPE_DEF
#define SX_TYPE_DEF

#include "config.h"

#if USE_LONG_DOUBLE
typedef long double              SX_FLT;
#else
typedef double                   SX_FLT;
#endif

#if USE_LONG_LONG
typedef signed long long int     SX_INT;
#elif USE_LONG
typedef signed long int          SX_INT;
#else
typedef signed int               SX_INT;
#endif

/**
 * \struct SX_MAT
 * \brief Sparse matrix of SX_FLT type in CSR format
 *
 * CSR Format (Ap, Aj, Ax) in SX_FLT
 *
 * \note The starting index of A is 0.
 */
typedef struct SX_MAT_
{
    SX_INT num_rows;
    SX_INT num_cols;
    SX_INT num_nnzs;

    SX_INT *Ap;
    SX_INT *Aj;
    SX_FLT *Ax;

} SX_MAT;

/**
 * \struct SX_MAT
 * \brief Sparse matrix of SX_INT type in CSR format
 *
 * CSR Format (Ap, Aj, Ax) in integer
 *
 * \note The starting index of A is 0.
 */
typedef struct SX_IMAT_
{
    SX_INT num_rows;
    SX_INT num_cols;
    SX_INT num_nnzs;

    SX_INT *Ap;
    SX_INT *Aj;
    SX_INT *Ax;

} SX_IMAT;

/**
 * \struct SX_VEC
 * \brief Vector with n entries of SX_FLT type
 */
typedef struct SX_VEC_
{
    SX_INT n;
    SX_FLT *d;

} SX_VEC;

/**
 * \struct SX_IVEC
 * \brief Vector with n entries of SX_INT type
 */
typedef struct SX_IVEC_
{
    SX_INT n;
    SX_INT *d;

} SX_IVEC;

/**
 * \brief Definition of standard smoother types
 */
typedef enum SX_SM_TYPE_
{
    SX_SM_JACOBI    = 1,  /**< Jacobi smoother */
    SX_SM_GS        = 2,  /**< Gauss-Seidel smoother */
    SX_SM_SGS       = 3,  /**< Symmetric Gauss-Seidel smoother */
    SX_SM_SOR       = 4,  /**< SOR smoother */
    SX_SM_SSOR      = 5,  /**< SSOR smoother */
    SX_SM_GSOR      = 6,  /**< GS + SOR smoother */
    SX_SM_SGSOR     = 7,  /**< SGS + SSOR smoother */
    SX_SM_POLY      = 8,  /**< Polynomial smoother */
    SX_SM_L1DIAG    = 9,  /**< L1 norm diagonal scaling smoother */

} SX_SM_TYPE;

/**
 * \brief Definition of coarsening types
 */
typedef enum SX_COARSEN_TYPE_
{
    SX_COARSE_RS      = 1,  /**< Classical */
    SX_COARSE_RSP     = 2,  /**< Classical, with positive offdiags */

} SX_COARSEN_TYPE;

/**
 * \brief Definition of interpolation types
 */
typedef enum SX_INTERP_TYPE_
{
    SX_INTERP_DIR     = 1,  /**< Direct interpolation */
    SX_INTERP_STD     = 2,  /**< Standard interpolation */

} SX_INTERP_TYPE;

typedef struct SX_SMTR_
{
    SX_SM_TYPE smoother;

    SX_MAT *A;
    SX_VEC *b;
    SX_VEC *x;

    SX_FLT relax;
    SX_INT nsweeps;
    SX_INT istart;
    SX_INT iend;
    SX_INT istep;
    SX_INT ndeg;
    SX_INT cf_order;
    SX_INT *ordering;

} SX_SMTR;

/**
 * \brief Definition of coarsest solver types
 */
typedef enum SX_CSSOLVE_TYPE_
{
    SX_CSSOLVE_LU        = 1,  /**< LU */
    SX_CSSOLVE_GMRES     = 2,  /**< gmres */
    SX_CSSOLVE_CG        = 3,  /**< cg */
    SX_CSSOLVE_BICGSTAB  = 4,  /**< bicgstab */

} SX_CSSOLVE_TYPE;

/**
 * \struct SX_AMG_PARS
 * \brief Parameters for AMG solver
 *
 */
typedef struct SX_AMG_PARS_
{
    SX_INT verb;

    SX_INT cycle_itr;            /** type of AMG cycle, 0, 1 is for V, others for W */
    SX_FLT tol;                  /** stopping tolerance for AMG solver */
    SX_FLT ctol;                 /** stopping tolerance for coarsest solver */
    SX_INT maxit;                /** max number of iterations of AMG */

    SX_COARSEN_TYPE cs_type;     /** coarsening type */
    SX_INT max_levels;           /** max number of levels of AMG */
    SX_INT coarse_dof;           /** max number of coarsest level DOF */
    SX_CSSOLVE_TYPE cssolve_type;

    SX_SM_TYPE smoother;         /** smoother type */
    SX_FLT relax;                /** relax parseter for SOR smoother */
    SX_INT cf_order;             /** False: nature order True: C/F order */
    SX_INT pre_iter;             /** number of presmoothers */
    SX_INT post_iter;            /** number of postsmoothers */
    SX_INT poly_deg;             /** degree of the polynomial smoother */

    SX_INTERP_TYPE interp_type;  /** interpolation type */
    SX_FLT strong_threshold;     /** strong connection threshold for coarsening */
    SX_FLT max_row_sum;          /** maximal row sum parseter */
    SX_FLT trunc_threshold;      /** truncation threshold */

} SX_AMG_PARS;

/* return values */
typedef struct SX_RTN_
{
    SX_FLT ares;     /* absolute residual */
    SX_FLT rres;     /* relative residual */
    SX_INT nits;     /* number of iterations */

} SX_RTN;

typedef struct SX_AMG_COMP_
{
    SX_MAT A;                    /** pointer to the matrix at level level_num */
    SX_MAT R;                    /** restriction operator at level level_num */
    SX_MAT P;                    /** prolongation operator at level level_num */
    SX_VEC b;                    /** pointer to the right-hand side at level level_num */
    SX_VEC x;                    /** pointer to the iterative solution at level level_num */
    SX_IVEC cfmark;              /** pointer to the CF marker at level level_num */

    SX_VEC wp;                   /** cache work space */

    SX_FLT *LU;
    int *pvt;

} SX_AMG_COMP;

/** 
 * \struct SX_AMG 
 * \brief Data for AMG solvers
 * 
 */ 
typedef struct SX_AMG
{
    SX_INT num_levels;           /** number of levels */
    SX_AMG_COMP *cg;             /* coarser grids */

    SX_AMG_PARS pars;            /* AMG parameters */

    SX_RTN rtn;                  /* return values */

} SX_AMG;

typedef struct SX_KRYLOV_
{
    SX_FLT tol;
    SX_MAT *A;
    SX_VEC *b;
    SX_VEC *u;
    SX_INT restart;
    SX_INT maxit;
    SX_INT verb;

    SX_RTN rtn;

} SX_KRYLOV;

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

#endif
