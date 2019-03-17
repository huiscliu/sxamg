
#include "amg-utils.h"
#include "internal.h"
#include <assert.h>

/**
 * \fn SX_AMG * sx_amg_data_create(SX_AMG_PARS *pars)
 *
 * \brief Create and initialize SX_AMG for classical AMG
 *
 * \pars pars AMG parameters
 *
 * \return Pointer to the SX_AMG data structure
 *
 */
SX_AMG sx_amg_data_create(SX_AMG_PARS *pars)
{
    SX_AMG mg;

    assert(pars != NULL);
    assert(pars->max_levels > 0);

    bzero(&mg, sizeof(mg));
    mg.cg = sx_calloc(pars->max_levels, sizeof(*mg.cg));

    /* record pars */
    mg.pars = *pars;

    return mg;
}

/**
 * \fn void sx_amg_data_destroy(SX_AMG **mgp)
 *
 * \brief Free SX_AMG data memeory space
 *
 * \pars mgp Pointer to pointer of SX_AMG
 *
 */
void sx_amg_data_destroy(SX_AMG *mg)
{
    SX_INT max_levels;
    SX_INT i;

    if (mg == NULL) return;

    max_levels = SX_MAX(1, mg->num_levels);

    for (i = 0; i < max_levels; ++i) {
        sx_mat_destroy(&mg->cg[i].A);
        sx_mat_destroy(&mg->cg[i].P);
        sx_mat_destroy(&mg->cg[i].R);

        if (i > 0) {
            sx_vec_destroy(&mg->cg[i].b);
            sx_vec_destroy(&mg->cg[i].x);
        }

        sx_vec_destroy(&mg->cg[i].wp);
        sx_ivec_destroy(&mg->cg[i].cfmark);
    }

    sx_free(mg->cg);
    bzero(mg, sizeof(*mg));
}

/**
 * \fn void void sx_amg_complexity_print (SX_AMG *mg)
 *
 * \brief Print complexities of AMG method
 *
 * \pars mg Multilevel hierachy for AMG
 *
 */
void sx_amg_complexity_print(SX_AMG *mg)
{
    const SX_INT max_levels = mg->num_levels;
    SX_INT level;
    SX_FLT gridcom = 0.0, opcom = 0.0;

    sx_printf("-----------------------------------------------------------\n");
    sx_printf("  Level   Num of rows   Num of nonzeros   Avg. NNZ / row   \n");
    sx_printf("-----------------------------------------------------------\n");

    for (level = 0; level < max_levels; ++level) {
        SX_FLT AvgNNZ = (SX_FLT) mg->cg[level].A.num_nnzs / mg->cg[level].A.num_rows;

        sx_printf("%5"dFMT" %13"dFMT" %17"dFMT" %14.2"fFMTf"\n", level, mg->cg[level].A.num_rows,
                mg->cg[level].A.num_nnzs, AvgNNZ);

        gridcom += mg->cg[level].A.num_rows;
        opcom += mg->cg[level].A.num_nnzs;
    }

    sx_printf("-----------------------------------------------------------\n");

    gridcom /= mg->cg[0].A.num_rows;
    opcom /= mg->cg[0].A.num_nnzs;
    sx_printf("  Grid complexity = %.3"fFMTf"  |", gridcom);
    sx_printf("  Operator complexity = %.3"fFMTf"\n", opcom);

    sx_printf("-----------------------------------------------------------\n");
}

/**
 * \fn void sx_amg_pars_init(SX_AMG_PARS *pars)
 *
 * \brief Initialize AMG parseters
 *
 * \pars pars Parameters for AMG
 *
 */
void sx_amg_pars_init(SX_AMG_PARS *pars)
{
    pars->verb = 0;

    pars->smoother = SX_SM_GS;

    pars->maxit = 100;
    pars->tol = 1e-6;
    pars->ctol = 1e-7;

    pars->max_levels = 30;
    pars->coarse_dof = MIN_CDOF;

    /* default v cycle */
    pars->cycle_itr = 1;

    pars->cf_order = 1;
    pars->pre_iter = 2;
    pars->post_iter = 2;
    pars->relax = 1.0;
    pars->poly_deg = 3;

    pars->cs_type = SX_COARSE_RSP;

    pars->interp_type = SX_INTERP_STD;
    pars->max_row_sum = 0.9;
    pars->strong_threshold = 0.3;
    pars->trunc_threshold = 0.2;
}

/**
 * \fn void sx_amg_pars_print(SX_AMG_PARS *pars)
 *
 * \brief Print out AMG parseters
 *
 * \pars pars Parameters for AMG
 *
 */
void sx_amg_pars_print(SX_AMG_PARS *pars)
{
    assert(pars != NULL);

    sx_printf("\n               AMG Parameters \n");
    sx_printf("-----------------------------------------------------------\n");

    sx_printf("AMG print level:                   %"dFMT"\n", pars->verb);
    sx_printf("AMG max num of iter:               %"dFMT"\n", pars->maxit);
    sx_printf("AMG tol:                           %"fFMTg"\n", pars->tol);
    sx_printf("AMG ctol:                          %"fFMTg"\n", pars->ctol);
    sx_printf("AMG max levels:                    %"dFMT"\n", pars->max_levels);
    sx_printf("AMG cycle type:                    %"dFMT"\n", pars->cycle_itr);
    sx_printf("AMG smoother type:                 %"dFMT"\n", pars->smoother);
    sx_printf("AMG smoother order:                %"dFMT"\n", pars->cf_order);
    sx_printf("AMG num of presmoothing:           %"dFMT"\n", pars->pre_iter);
    sx_printf("AMG num of postsmoothing:          %"dFMT"\n", pars->post_iter);

    switch(pars->smoother) {
        case SX_SM_SOR:
        case SX_SM_SSOR:
        case SX_SM_GSOR:
        case SX_SM_SGSOR:
            sx_printf("AMG relax factor:                  %.4"fFMTf"\n", pars->relax);
            break;

        case SX_SM_POLY:
            sx_printf("AMG polynomial smoother degree:    %"dFMT"\n", pars->poly_deg);
            break;

        default:
            break;
    }

    sx_printf("AMG coarsening type:               %"dFMT"\n", pars->cs_type);
    sx_printf("AMG interpolation type:            %"dFMT"\n", pars->interp_type);
    sx_printf("AMG dof on coarsest grid:          %"dFMT"\n", pars->coarse_dof);
    sx_printf("AMG strong threshold:              %.4"fFMTf"\n", pars->strong_threshold);
    sx_printf("AMG truncation threshold:          %.4"fFMTf"\n", pars->trunc_threshold);
    sx_printf("AMG max row sum:                   %.4"fFMTf"\n", pars->max_row_sum);

    sx_printf("-----------------------------------------------------------\n");
}

/**
 * \fn void sx_print_itinfo (const SX_INT verb, const SX_INT stop_type, const SX_INT iter,
 *                        const SX_FLT relres, const SX_FLT absres, const SX_FLT factor)
 *
 * \brief Print out iteration information for iterative solvers
 *
 * \pars verb     Level for output
 * \pars stop_type  Type of stopping criteria
 * \pars iter       Number of iterations
 * \pars relres     Relative residual of different kinds
 * \pars absres     Absolute residual of different kinds
 * \pars factor     Contraction factor
 *
 */
void sx_print_itinfo(const SX_INT verb, const SX_INT stop_type, const SX_INT iter, const SX_FLT relres,
        const SX_FLT absres, const SX_FLT factor)
{
    if (verb >= 2) {

        if (iter > 0) {
            sx_printf("%6"dFMT" | %13.6"fFMTe"   | %13.6"fFMTe"  | %10.4"fFMTf"\n", iter, relres, absres,
                    factor);
        }
        else {                  // iter = 0: initial guess
            sx_printf("-----------------------------------------------------------\n");

            switch (stop_type) {
                case STOP_REL_RES:
                    sx_printf("It Num |   ||r||/||b||   |     ||r||      |  Conv. Factor\n");
                    break;

                case STOP_REL_PRECRES:
                    sx_printf("It Num | ||r||_B/||b||_B |    ||r||_B     |  Conv. Factor\n");
                    break;

                case STOP_MOD_REL_RES:
                    sx_printf("It Num |   ||r||/||x||   |     ||r||      |  Conv. Factor\n");
                    break;
            }

            sx_printf("-----------------------------------------------------------\n");
            sx_printf("%6"dFMT" | %13.6"fFMTe"   | %13.6"fFMTe"  |     -.-- \n", iter, relres, absres);
        }
    }
}
