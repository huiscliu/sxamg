
#include "amg-cycle.h"
#include "internal.h"

/**
 * \fn void sx_amg_cycle(SX_AMG *mg, SX_AMG_PARS *pars)
 *
 * \brief Solve Ax=b with non-recursive multigrid cycle
 *
 * \pars mg    Pointer to AMG data: SX_AMG
 * \pars pars  Pointer to AMG parseters: SX_AMG_PARS
 *
 */
void sx_amg_cycle(SX_AMG *mg)
{
    SX_INT verb = mg->pars.verb;
    SX_INT cycle_itr = mg->pars.cycle_itr;
    SX_INT nl = mg->num_levels;
    SX_FLT tol = mg->pars.ctol;

    SX_FLT alpha = 1.0;
    SX_INT num_lvl[MAX_AMG_LVL] = {0}, l = 0;

    if (tol > mg->pars.tol)  tol = mg->pars.tol * 0.1;
    if (cycle_itr <= 0) cycle_itr = 1;

ForwardSweep:
    while (l < nl - 1) {
        SX_SMTR s;

        num_lvl[l]++;

        // pre-smoothing
        s.smoother = mg->pars.smoother;
        s.A = &mg->cg[l].A;
        s.b = &mg->cg[l].b;
        s.x = &mg->cg[l].x;
        s.nsweeps = mg->pars.pre_iter;
        s.istart = 0;
        s.iend = mg->cg[l].A.num_rows - 1;
        s.istep = 1;
        s.relax = mg->pars.relax;
        s.ndeg = mg->pars.poly_deg;
        s.cf_order = mg->pars.cf_order;
        s.ordering = mg->cg[l].cfmark.d;

        sx_amg_smoother_pre(&s);

        // form residual r = b - A x
        sx_blas_array_cp(mg->cg[l].A.num_rows, mg->cg[l].b.d, mg->cg[l].wp.d);
        sx_blas_mv_amxpy(-1.0, &mg->cg[l].A, &mg->cg[l].x, &mg->cg[l].wp);

        // restriction r1 = R*r0
        sx_blas_mv_mxy(&mg->cg[l].R, &mg->cg[l].wp, &mg->cg[l + 1].b);

        // prepare for the next level
        l++;
        sx_vec_set_value(&mg->cg[l].x, 0.0);
    }

    // call the coarse space solver:
    sx_amg_coarest_solve(&mg->cg[nl - 1].A, &mg->cg[nl - 1].b, &mg->cg[nl - 1].x, tol, verb);

    // BackwardSweep:
    while (l > 0) {
        SX_SMTR s;

        l--;

        // prolongation u = u + alpha*P*e1
        sx_blas_mv_amxpy(alpha, &mg->cg[l].P, &mg->cg[l + 1].x, &mg->cg[l].x);

        // post-smoothing
        s.smoother = mg->pars.smoother;
        s.A = &mg->cg[l].A;
        s.b = &mg->cg[l].b;
        s.x = &mg->cg[l].x;
        s.nsweeps = mg->pars.post_iter;
        s.istart = 0;
        s.iend = mg->cg[l].A.num_rows - 1;
        s.istep = -1;
        s.relax = mg->pars.relax;
        s.ndeg = mg->pars.poly_deg;
        s.cf_order = mg->pars.cf_order;
        s.ordering = mg->cg[l].cfmark.d;

        sx_amg_smoother_post(&s);

        if (num_lvl[l] < cycle_itr)
            break;
        else
            num_lvl[l] = 0;
    }

    if (l > 0) goto ForwardSweep;
}
