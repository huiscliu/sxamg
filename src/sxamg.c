
#include "sxamg.h"
#include "internal.h"

SX_RTN sx_solver_amg_solve(SX_AMG *mg, SX_VEC *x, SX_VEC *b)
{
    SX_MAT *ptrA;
    SX_VEC *r;

    SX_INT verb;
    SX_INT MaxIt;
    SX_FLT tol;
    SX_FLT sumb; // L2norm(b)

    SX_FLT solve_start, solve_end;
    SX_FLT relres1 = 1., absres0, absres, factor;
    SX_INT iter = 0;
    SX_RTN rtn;

    /* check pars */
    assert(mg != NULL);
    assert(x != NULL);
    assert(b != NULL);

    ptrA = &mg->cg[0].A;
    r = &mg->cg[0].wp;

    verb = mg->pars.verb;
    MaxIt = mg->pars.maxit;
    tol = mg->pars.tol;
    sumb = sx_blas_vec_norm2(b); // L2norm(b)
    absres0 = sumb;

    solve_start = sx_get_time();

    // Print iteration information if needed
    sx_print_itinfo(verb, STOP_REL_RES, iter, 1.0, sumb, 0.0);

    /* init, to make compiler happy */
    rtn.ares = 0;
    rtn.rres = 0;
    rtn.nits = 0;

    if (fabs(sumb) == 0.) {
        sx_vec_set_value(x, 0);
        mg->rtn = rtn;

        return rtn;
    }

    /* set x and b */
    mg->cg[0].x = *x;
    mg->cg[0].b = *b;

    // MG solver here
    while ((++iter <= MaxIt)) {

        sx_amg_cycle(mg);

        // Form residual r = b - A*x
        sx_vec_cp(b, r);
        sx_blas_mv_amxpy(-1.0, ptrA, x, r);

        // Compute norms of r and convergence factor
        absres = sx_blas_vec_norm2(r);  // residual ||r||
        relres1 = absres / sumb;        // relative residual ||r||/||b||
        factor = absres / absres0;      // contraction factor
        absres0 = absres;               // prepare for next iteration

        // Print iteration information if needed
        sx_print_itinfo(verb, STOP_REL_RES, iter, relres1, absres, factor);

        /* save convergence info */
        rtn.ares = absres;
        rtn.rres = relres1;
        rtn.nits = iter;
        mg->rtn = rtn;

        // Check convergence
        if (relres1 < tol) break;
    }

    if (verb > 0) {
        solve_end = sx_get_time();
        sx_printf("AMG solve time: %"fFMTg" s\n", solve_end - solve_start);
    }

    return rtn;
}

/**
 * \fn void sx_solver_amg (SX_MAT *A, SX_VEC *b, SX_VEC *x, SX_AMG_PARS *pars)
 *
 * \brief Solve Ax = b by algebraic multigrid methods
 *
 * \pars A     Pointer to SX_MAT: the coefficient matrix
 * \pars x     Pointer to SX_VEC: the unknowns
 * \pars b     Pointer to SX_VEC: the right hand side
 * \pars pars  Pointer to SX_AMG_PARS: AMG parseters
 *
 * \note Refer to "Multigrid"
 *       by U. Trottenberg, C. W. Oosterlee and A. Schuller
 *       Appendix A.7 (by A. Brandt, P. Oswald and K. Stuben)
 *       Academic Press Inc., San Diego, CA, 2001.
 *
 */
SX_RTN sx_solver_amg(SX_MAT *A, SX_VEC *x, SX_VEC *b, SX_AMG_PARS *pars)
{
    SX_INT verb;
    SX_INT nnz, m, n;
    SX_RTN rtn;
    SX_AMG_PARS npars;

    SX_AMG mg;
    SX_FLT AMG_start, AMG_end;
    SX_FLT sumb;

    assert(A != NULL);
    assert(b != NULL);
    assert(x != NULL);

    sumb = sx_blas_vec_norm2(b);

    if (pars == NULL) {
        sx_amg_pars_init(&npars);
        pars = &npars;
    }

    verb = pars->verb;
    if (fabs(sumb) == 0.) {
        sx_vec_set_value(x, 0);
        rtn.ares = 0;
        rtn.rres = 0;
        rtn.nits = 0;

        sx_print_itinfo(verb, STOP_REL_RES, 0, 0., sumb, 0.0);

        return rtn;
    }

    nnz = A->num_nnzs;
    m = A->num_rows;
    n = A->num_cols;

    if (verb > 0) AMG_start = sx_get_time();

    // check matrix data
    if (m != n) {
        sx_printf("### ERROR: A is not a square matrix!\n");
        sx_exit_on_errcode(ERROR_MAT_SIZE, __FUNCTION__);
    }

    if (nnz <= 0) {
        sx_printf("### ERROR: A has no nonzero entries!\n");
        sx_exit_on_errcode(ERROR_MAT_SIZE, __FUNCTION__);
    }

    // Step 1: AMG setup phase
    sx_amg_setup(&mg, A, pars);

    // Step 2: AMG solve phase
    rtn = sx_solver_amg_solve(&mg, x, b);

    // clean-up memory
    sx_amg_data_destroy(&mg);

    // print out CPU time if needed
    if (verb > 0) {
        AMG_end = sx_get_time();
        sx_printf("AMG totally time: %"fFMTg" s\n", AMG_end - AMG_start);
    }

    return rtn;
}
