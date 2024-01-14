
#include "amg-coarse-solver.h"
#include "amg-krylov.h"
#include "internal.h"

/**
 * \fn static void sx_amg_coarest_solve(SX_AMG_COMP *cg, const SX_FLT ctol, const SX_INT verb)
 *
 * \brief Iterative on the coarset level
 *
 * \pars  cg     pointer to matrix data
 * \pars  ctol   tolerance for the coarsest level
 * \pars  cssolve_type solver type for coarsest level
 * \pars  verb   level of output
 *
 */
void sx_amg_coarest_solve(SX_AMG_COMP *cg, const SX_FLT ctol, SX_CSSOLVE_TYPE cssolve_type, const SX_INT verb)
{
    SX_MAT *A = &cg->A;
    SX_VEC *b = &cg->b;
    SX_VEC *x = &cg->x;
    const SX_INT n = A->num_rows;
    const SX_INT maxit = SX_MAX(250, SX_MIN(n * n, 1000));
    SX_RTN rtn ;
    SX_KRYLOV ks;

    switch(cssolve_type) {
        case SX_CSSOLVE_CG:
            /* cg */
            ks.A = A;
            ks.b = b;
            ks.u = x;
            ks.tol = ctol;
            ks.maxit = maxit;
            ks.verb = 0;

            rtn = sx_solver_cg_itnl(&ks, NULL);

            if (rtn.nits >= maxit && verb >= 2) {
                sx_printf("### WARNING: Coarse level solver failed to converge!\n");
            }

            break;

        case SX_CSSOLVE_GMRES:
            /* cg */
            ks.A = A;
            ks.b = b;
            ks.u = x;
            ks.tol = ctol;
            ks.maxit = maxit;
            ks.verb = 0;

            /* try GMRES if cg fails */
            ks.restart = MAX_RESTART;
            rtn = sx_solver_gmres_itnl(&ks, NULL);

            if (rtn.nits >= maxit && verb >= 2) {
                sx_printf("### WARNING: Coarse level solver failed to converge!\n");
            }

            break;

        case SX_CSSOLVE_BICGSTAB:
            /* cg */
            ks.A = A;
            ks.b = b;
            ks.u = x;
            ks.tol = ctol;
            ks.maxit = maxit;
            ks.verb = 0;

            rtn = sx_solver_bicgstab_itnl(&ks, NULL);

            if (rtn.nits >= maxit && verb >= 2) {
                sx_printf("### WARNING: Coarse level solver failed to converge!\n");
            }

            break;

        case SX_CSSOLVE_LU:

            /* init */
            if (cg->LU == NULL) {
                SX_INT i, k, j, ibegin, iend;

                cg->LU = sx_calloc(n * (size_t) n, sizeof(*cg->LU));
                cg->pvt = sx_calloc(n, sizeof(*cg->pvt));

                /* init LU */
                for (i = 0; i < n; ++i) {
                    ibegin = A->Ap[i];
                    iend = A->Ap[i + 1];
                    for (k = ibegin; k < iend; ++k) {
                        j = A->Aj[k];

                        cg->LU[i * n + j] = A->Ax[k];
                    }
                }

                /* LU factorization */
                sx_mat_dense_lu(n, cg->LU, cg->pvt);
            }

            /* init x */
            memcpy(x->d, b->d, n * sizeof(*x->d));

            /* solve */
            sx_mat_dense_sv(n, cg->LU, cg->pvt, 1, x->d);

            break;
    }
}
