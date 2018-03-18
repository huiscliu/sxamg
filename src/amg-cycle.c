
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
    SX_INT nl = mg[0].num_levels;
    SX_FLOAT tol = mg->pars.ctol;

    SX_FLOAT alpha = 1.0;
    SX_INT num_lvl[MAX_AMG_LVL] = {0}, l = 0;

    if (tol > mg->pars.tol)  tol = mg->pars.tol * 0.1;
    if (cycle_itr <= 0) cycle_itr = 1;

ForwardSweep:
    while (l < nl - 1) {
        SX_SMTR s;

        num_lvl[l]++;

        // pre-smoothing
        s.smoother = mg->pars.smoother;
        s.A = &mg[l].A;
        s.b = &mg[l].b;
        s.x = &mg[l].x;
        s.nsweeps = mg->pars.pre_iter;
        s.istart = 0;
        s.iend = mg[l].A.num_rows - 1;
        s.istep = 1;
        s.relax = mg->pars.relax;
        s.ndeg = mg->pars.poly_deg;
        s.cf_order = mg->pars.cf_order;
        s.ordering = mg[l].cfmark.d;

        sx_amg_smoother_pre(&s);

        // form residual r = b - A x
        sx_blas_array_cp(mg[l].A.num_rows, mg[l].b.d, mg[l].wp.d);
        sx_blas_mat_amxpy(-1.0, &mg[l].A, &mg[l].x, &mg[l].wp);

        // restriction r1 = R*r0
        sx_blas_mat_mxy(&mg[l].R, &mg[l].wp, &mg[l + 1].b);

        // prepare for the next level
        l++;
        sx_vec_set_value(&mg[l].x, 0.0);
    }

    // call the coarse space solver:
    sx_amg_coarest_solve(&mg[nl - 1].A, &mg[nl - 1].b, &mg[nl - 1].x, tol, verb);

    // BackwardSweep:
    while (l > 0) {
        SX_SMTR s;

        l--;

        // prolongation u = u + alpha*P*e1
        sx_blas_mat_amxpy(alpha, &mg[l].P, &mg[l + 1].x, &mg[l].x);

        // post-smoothing
        s.smoother = mg->pars.smoother;
        s.A = &mg[l].A;
        s.b = &mg[l].b;
        s.x = &mg[l].x;
        s.nsweeps = mg->pars.post_iter;
        s.istart = 0;
        s.iend = mg[l].A.num_rows - 1;
        s.istep = -1;
        s.relax = mg->pars.relax;
        s.ndeg = mg->pars.poly_deg;
        s.cf_order = mg->pars.cf_order;
        s.ordering = mg[l].cfmark.d;

        sx_amg_smoother_post(&s);

        if (num_lvl[l] < cycle_itr)
            break;
        else
            num_lvl[l] = 0;
    }

    if (l > 0) goto ForwardSweep;
}
