
#include "amg-coarse-solver.h"
#include "internal.h"

static SX_INT sx_solver_cg(SX_KRYLOV *ks)
{
    SX_MAT *A = ks->A;
    SX_VEC *b = ks->b;
    SX_VEC *u = ks->u;
    SX_FLT tol = ks->tol;
    SX_INT maxit = ks->maxit;
    SX_INT stop_type = ks->stop_type;
    SX_INT verb = ks->verb;

    SX_INT MaxStag = MAX_STAG, MaxRestartStep = MAX_RESTART;
    SX_INT m = b->n;
    SX_FLT maxdiff = tol * 1e-4;       // staganation tolerance
    SX_FLT sol_inf_tol = SMALLFLOAT;   // infinity norm tolerance
    SX_INT iter = 0, stag = 1, more_step = 1, restart_step = 1;
    SX_FLT absres0 = BIGFLOAT, absres = BIGFLOAT;
    SX_FLT relres = BIGFLOAT, normu = BIGFLOAT, normr0 = BIGFLOAT;
    SX_FLT reldiff, factor, infnormu;
    SX_FLT alpha, beta, temp1, temp2;

    SX_INT iter_best = 0;              // initial best known iteration
    SX_FLT absres_best = BIGFLOAT;   // initial best known residual

    // allocate temp memory (need 5*m SX_FLT numbers)
    SX_FLT *work = (SX_FLT *) sx_calloc(5 * m, sizeof(SX_FLT));
    SX_FLT *p = work, *z = work + m, *r = z + m, *t = r + m, *u_best = t + m;
    SX_VEC vr;

    vr.n = b->n;
    vr.d = r;

    // r = b-A*u
    sx_blas_array_cp(m, b->d, r);
    sx_blas_mv_amxpy(-1.0, A, u, &vr);
    sx_blas_array_cp(m, r, z);

    // compute initial residuals
    switch (stop_type) {
        case STOP_REL_RES:
            absres0 = sx_blas_array_norm2(m, r);
            normr0 = SX_MAX(SMALLFLOAT, absres0);
            relres = absres0 / normr0;
            break;

        case STOP_REL_PRECRES:
            absres0 = sqrt(sx_blas_array_dot(m, r, z));
            normr0 = SX_MAX(SMALLFLOAT, absres0);
            relres = absres0 / normr0;
            break;

        case STOP_MOD_REL_RES:
            absres0 = sx_blas_array_norm2(m, r);
            normu = SX_MAX(SMALLFLOAT, sx_blas_array_norm2(m, u->d));
            relres = absres0 / normu;
            break;

        default:
            sx_printf("### ERROR: Unrecognized stopping type for %s!\n", __FUNCTION__);
            goto eofc;
    }

    // if initial residual is small, no need to iterate!
    if (relres < tol) goto eofc;

    // output iteration information if needed
    sx_print_itinfo(verb, stop_type, iter, relres, absres0, 0.0);
    sx_blas_array_cp(m, z, p);
    temp1 = sx_blas_array_dot(m, z, r);

    // main CG loop
    while (iter++ < maxit) {
        SX_VEC vp, vt;

        // t=A*p
        vp.n = b->n;
        vp.d = p;
        vt.n = b->n;
        vt.d = t;
        sx_blas_mv_mxy(A, &vp, &vt);

        // alpha_k=(z_{k-1},r_{k-1})/(A*p_{k-1},p_{k-1})
        temp2 = sx_blas_array_dot(m, t, p);
        if (SX_ABS(temp2) > SMALLFLOAT2) {
            alpha = temp1 / temp2;
        }
        else {                  // Possible breakdown
            goto RESTORE_BESTSOL;
        }

        // u_k=u_{k-1} + alpha_k*p_{k-1}
        sx_blas_array_axpy(m, alpha, p, u->d);

        // r_k=r_{k-1} - alpha_k*A*p_{k-1}
        sx_blas_array_axpy(m, -alpha, t, r);

        // compute residuals
        switch (stop_type) {
            case STOP_REL_RES:
                absres = sx_blas_array_norm2(m, r);
                relres = absres / normr0;
                break;
            case STOP_REL_PRECRES:
                // z = B(r)
                sx_blas_array_cp(m, r, z); /* No preconditioner */
                absres = sqrt(SX_ABS(sx_blas_array_dot(m, z, r)));
                relres = absres / normr0;
                break;
            case STOP_MOD_REL_RES:
                absres = sx_blas_array_norm2(m, r);
                relres = absres / normu;
                break;
        }

        // compute reducation factor of residual ||r||
        factor = absres / absres0;

        // output iteration information if needed
        sx_print_itinfo(verb, stop_type, iter, relres, absres, factor);

        // safety net check: save the best-so-far solution
        if (absres < absres_best - maxdiff) {
            absres_best = absres;
            iter_best = iter;
            sx_blas_array_cp(m, u->d, u_best);
        }

        // Check I: if soultion is close to zero, return ERROR_SOLVER_SOLSTAG
        infnormu = sx_blas_array_norminf(m, u->d);
        if (infnormu <= sol_inf_tol) {
            iter = ERROR_SOLVER_SOLSTAG;
            break;
        }

        //  Check II: if staggenated, try to restart
        normu = sx_blas_vec_norm2(u);

        // compute relative difference
        reldiff = SX_ABS(alpha) * sx_blas_array_norm2(m, p) / normu;
        if ((stag <= MaxStag) & (reldiff < maxdiff)) {

            sx_blas_array_cp(m, b->d, r);
            sx_blas_mv_amxpy(-1.0, A, u, &vr);

            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = sx_blas_array_norm2(m, r);
                    relres = absres / normr0;
                    break;

                case STOP_REL_PRECRES:
                    // z = B(r)
                    sx_blas_array_cp(m, r, z);     /* No preconditioner */
                    absres = sqrt(SX_ABS(sx_blas_array_dot(m, z, r)));
                    relres = absres / normr0;
                    break;

                case STOP_MOD_REL_RES:
                    absres = sx_blas_array_norm2(m, r);
                    relres = absres / normu;
                    break;
            }

            if (relres < tol) {
                break;
            }
            else {
                if (stag >= MaxStag) {
                    iter = ERROR_SOLVER_STAG;
                    break;
                }
                sx_blas_array_set(m, p, 0.0);
                ++stag;
                ++restart_step;
            }
        }                       // end of staggnation check!

        // Check III: prevent false convergence
        if (relres < tol) {
            // compute residual r = b - Ax again
            sx_blas_array_cp(m, b->d, r);
            sx_blas_mv_amxpy(-1.0, A, u, &vr);

            // compute residuals
            switch (stop_type) {
                case STOP_REL_RES:
                    absres = sx_blas_array_norm2(m, r);
                    relres = absres / normr0;
                    break;

                case STOP_REL_PRECRES:
                    // z = B(r)
                    sx_blas_array_cp(m, r, z);     /* No preconditioner */
                    absres = sqrt(SX_ABS(sx_blas_array_dot(m, z, r)));
                    relres = absres / normr0;
                    break;

                case STOP_MOD_REL_RES:
                    absres = sx_blas_array_norm2(m, r);
                    relres = absres / normu;
                    break;
            }

            // check convergence
            if (relres < tol) break;

            if (more_step >= MaxRestartStep) {
                iter = ERROR_SOLVER_TOLSMALL;
                break;
            }

            // prepare for restarting the method
            sx_blas_array_set(m, p, 0.0);
            ++more_step;
            ++restart_step;
        }

        // save residual for next iteration
        absres0 = absres;

        // compute z_k = B(r_k)
        if (stop_type != STOP_REL_PRECRES) {
            sx_blas_array_cp(m, r, z);     /* No preconditioner, B=I */
        }

        // compute beta_k = (z_k, r_k)/(z_{k-1}, r_{k-1})
        temp2 = sx_blas_array_dot(m, z, r);
        beta = temp2 / temp1;
        temp1 = temp2;

        // compute p_k = z_k + beta_k*p_{k-1}
        sx_blas_array_axpby(m, 1.0, z, beta, p);

    }

RESTORE_BESTSOL:
    if (iter != iter_best) {
        SX_VEC vb;

        vb.n = b->n;
        vb.d = u_best;

        // compute best residual
        sx_blas_array_cp(m, b->d, r);
        sx_blas_mv_amxpy(-1.0, A, &vb, &vr);

        switch (stop_type) {
            case STOP_REL_RES:
                absres_best = sx_blas_array_norm2(m, r);
                break;

            case STOP_REL_PRECRES:
                // z = B(r)
                sx_blas_array_cp(m, r, z); /* No preconditioner */
                absres_best = sqrt(SX_ABS(sx_blas_array_dot(m, z, r)));
                break;

            case STOP_MOD_REL_RES:
                absres_best = sx_blas_array_norm2(m, r);
                break;
        }

        if (absres > absres_best + maxdiff) {
            sx_blas_array_cp(m, u_best, u->d);
            relres = absres_best / normr0;
        }
    }

eofc:

    // clean up temp memory
    sx_free(work);

    if (iter > maxit)
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

static SX_INT sx_solver_gmres(SX_KRYLOV *ks)
{
    SX_MAT *A = ks->A;
    SX_VEC *b = ks->b;
    SX_VEC *x = ks->u;
    SX_FLT tol = ks->tol;
    SX_INT maxit = ks->maxit;
    SX_INT restart = ks->restart;
    SX_INT stop_type = ks->stop_type;
    SX_INT verb = ks->verb;

    const SX_INT n = b->n;
    const SX_INT MIN_ITER = 0;
    const SX_FLT maxdiff = tol * 1e-4;     // staganation tolerance
    const SX_FLT epsmac = SMALLFLOAT;
    SX_INT iter = 0;
    SX_INT restart1 = restart + 1;
    SX_INT i, j, k;

    SX_FLT r_norm, r_normb, gamma, t;
    SX_FLT normr0 = BIGFLOAT, absres = BIGFLOAT;
    SX_FLT relres = BIGFLOAT, normu = BIGFLOAT;

    SX_INT iter_best = 0;          // initial best known iteration
    SX_FLT absres_best = BIGFLOAT;       // initial best known residual

    // allocate temp memory (need about (restart+4)*n FLOAT numbers)
    SX_FLT *c = NULL, *s = NULL, *rs = NULL;
    SX_FLT *norms = NULL, *r = NULL, *w = NULL;
    SX_FLT *work = NULL, *x_best = NULL;
    SX_FLT **p = NULL, **hh = NULL;
    SX_VEC vs, vr;

    /* allocate memory and setup temp work space */
    work = (SX_FLT *) sx_calloc((restart + 4) * (restart + n) + 1, sizeof(SX_FLT));

    /* check whether memory is enough for GMRES */
    while ((work == NULL) && (restart > 5)) {
        restart = restart - 5;
        work = (SX_FLT *) sx_calloc((restart + 4) * (restart + n) + 1, sizeof(SX_FLT));
        sx_printf("### WARNING: GMRES restart number set to %"dFMT"!\n", restart);
        restart1 = restart + 1;
    }

    if (work == NULL) {
        sx_printf("### ERROR: No enough memory for GMRES %s : %s: %d !\n", __FILE__,
                __FUNCTION__, __LINE__);

        exit(ERROR_ALLOC_MEM);
    }

    p = (SX_FLT **) sx_calloc(restart1, sizeof(SX_FLT *));
    hh = (SX_FLT **) sx_calloc(restart1, sizeof(SX_FLT *));
    norms = (SX_FLT *) sx_calloc(maxit + 1, sizeof(SX_FLT));

    r = work;
    w = r + n;
    rs = w + n;
    c = rs + restart1;
    x_best = c + restart;
    s = x_best + n;
    vr.n = vs.n = b->n;

    for (i = 0; i < restart1; i++) p[i] = s + restart + i * n;

    for (i = 0; i < restart1; i++) hh[i] = p[restart] + n + i * restart;

    // r = b-A*x
    sx_blas_array_cp(n, b->d, p[0]);

    vs.d = p[0];
    sx_blas_mv_amxpy(-1.0, A, x, &vs);
    r_norm = sx_blas_array_norm2(n, p[0]);

    // compute initial residuals
    switch (stop_type) {
        case STOP_REL_RES:
            normr0 = SX_MAX(SMALLFLOAT, r_norm);
            relres = r_norm / normr0;
            break;

        case STOP_REL_PRECRES:
            sx_blas_array_cp(n, p[0], r);
            r_normb = sqrt(sx_blas_array_dot(n, p[0], r));
            normr0 = SX_MAX(SMALLFLOAT, r_normb);
            relres = r_normb / normr0;
            break;

        case STOP_MOD_REL_RES:
            normu = SX_MAX(SMALLFLOAT, sx_blas_array_norm2(n, x->d));
            normr0 = r_norm;
            relres = normr0 / normu;
            break;

        default:
            sx_printf("### ERROR: Unrecognized stopping type for %s!\n", __FUNCTION__);
            goto eofc;
    }

    // if initial residual is small, no need to iterate!
    if (relres < tol) goto eofc;

    // output iteration information if needed
    sx_print_itinfo(verb, stop_type, 0, relres, normr0, 0.0);

    // store initial residual
    norms[0] = relres;

    /* outer iteration cycle */
    while (iter < maxit) {
        rs[0] = r_norm;
        t = 1.0 / r_norm;
        sx_blas_array_ax(n, t, p[0]);

        /* RESTART CYCLE (right-preconditioning) */
        i = 0;
        while (i < restart && iter < maxit) {
            SX_VEC vp;

            i++;
            iter++;

            vr.n = b->n;
            vp.n = b->n;
            vr.d = r;
            vp.d = p[i];

            sx_blas_array_cp(n, p[i - 1], r);
            sx_blas_mv_mxy(A, &vr, &vp);

            /* modified Gram_Schmidt */
            for (j = 0; j < i; j++) {
                hh[j][i - 1] = sx_blas_array_dot(n, p[j], p[i]);
                sx_blas_array_axpy(n, -hh[j][i - 1], p[j], p[i]);
            }
            t = sx_blas_array_norm2(n, p[i]);
            hh[i][i - 1] = t;
            if (t != 0.0) {
                t = 1.0 / t;
                sx_blas_array_ax(n, t, p[i]);
            }

            for (j = 1; j < i; ++j) {
                t = hh[j - 1][i - 1];
                hh[j - 1][i - 1] = s[j - 1] * hh[j][i - 1] + c[j - 1] * t;
                hh[j][i - 1] = -s[j - 1] * t + c[j - 1] * hh[j][i - 1];
            }

            t = hh[i][i - 1] * hh[i][i - 1];
            t += hh[i - 1][i - 1] * hh[i - 1][i - 1];

            gamma = sqrt(t);
            if (gamma == 0.0) gamma = epsmac;

            c[i - 1] = hh[i - 1][i - 1] / gamma;
            s[i - 1] = hh[i][i - 1] / gamma;
            rs[i] = -s[i - 1] * rs[i - 1];
            rs[i - 1] = c[i - 1] * rs[i - 1];
            hh[i - 1][i - 1] = s[i - 1] * hh[i][i - 1] + c[i - 1] * hh[i - 1][i - 1];

            absres = r_norm = fabs(rs[i]);

            relres = absres / normr0;

            norms[iter] = relres;

            // output iteration information if needed
            sx_print_itinfo(verb, stop_type, iter, relres, absres, norms[iter] / norms[iter - 1]);

            // should we exit the restart cycle
            if (relres <= tol && iter >= MIN_ITER) break;

        }

        /* compute solution, first solve upper triangular system */
        rs[i - 1] = rs[i - 1] / hh[i - 1][i - 1];
        for (k = i - 2; k >= 0; k--) {
            t = 0.0;
            for (j = k + 1; j < i; j++) t -= hh[k][j] * rs[j];

            t += rs[k];
            rs[k] = t / hh[k][k];
        }

        sx_blas_array_cp(n, p[i - 1], w);
        sx_blas_array_ax(n, rs[i - 1], w);

        for (j = i - 2; j >= 0; j--) sx_blas_array_axpy(n, rs[j], p[j], w);

        /* apply the preconditioner */
        sx_blas_array_cp(n, w, r);
        sx_blas_array_axpy(n, 1.0, r, x->d);

        if (absres < absres_best - maxdiff) {
            absres_best = absres;
            iter_best = iter;
            sx_blas_array_cp(n, x->d, x_best);
        }

        // Check: prevent false convergence
        if (relres <= tol && iter >= MIN_ITER) {

            sx_blas_array_cp(n, b->d, r);

            vs.d = r;
            sx_blas_mv_amxpy(-1.0, A, x, &vs);

            r_norm = sx_blas_array_norm2(n, r);

            switch (stop_type) {
                case STOP_REL_RES:
                    absres = r_norm;
                    relres = absres / normr0;
                    break;

                case STOP_REL_PRECRES:
                    sx_blas_array_cp(n, r, w);
                    absres = sqrt(sx_blas_array_dot(n, w, r));
                    relres = absres / normr0;
                    break;

                case STOP_MOD_REL_RES:
                    absres = r_norm;
                    normu = SX_MAX(SMALLFLOAT, sx_blas_array_norm2(n, x->d));
                    relres = absres / normu;
                    break;
            }

            norms[iter] = relres;

            if (relres <= tol) {
                break;
            }
            else {
                // Need to restart
                sx_blas_array_cp(n, r, p[0]);
                i = 0;
            }

        }

        /* compute residual vector and continue loop */
        for (j = i; j > 0; j--) {
            rs[j - 1] = -s[j - 1] * rs[j];
            rs[j] = c[j - 1] * rs[j];
        }

        if (i) sx_blas_array_axpy(n, rs[i] - 1.0, p[i], p[i]);

        for (j = i - 1; j > 0; j--)
            sx_blas_array_axpy(n, rs[j], p[j], p[i]);

        if (i) {
            sx_blas_array_axpy(n, rs[0] - 1.0, p[0], p[0]);
            sx_blas_array_axpy(n, 1.0, p[i], p[0]);
        }
    }                           /* end of main while loop */

    if (iter != iter_best) {
        // compute best residual
        sx_blas_array_cp(n, b->d, r);

        vs.d = x_best;
        vr.d = r;
        sx_blas_mv_amxpy(-1.0, A, &vs, &vr);

        switch (stop_type) {
            case STOP_REL_RES:
                absres_best = sx_blas_array_norm2(n, r);
                break;

            case STOP_REL_PRECRES:
                // z = B(r)
                sx_blas_array_cp(n, r, w); /* No preconditioner */
                absres_best = sqrt(SX_ABS(sx_blas_array_dot(n, w, r)));
                break;

            case STOP_MOD_REL_RES:
                absres_best = sx_blas_array_norm2(n, r);
                break;
        }

        if (absres > absres_best + maxdiff) {
            sx_blas_array_cp(n, x_best, x->d);
            relres = absres_best / normr0;
        }
    }

eofc:
    sx_free(work);
    sx_free(p);
    sx_free(hh);
    sx_free(norms);

    if (iter >= maxit)
        return ERROR_SOLVER_MAXIT;
    else
        return iter;
}

/**
 * \fn static void sx_amg_coarest_solve(SX_MAT *A, SX_VEC *b, SX_VEC *x,
 *                                      const SX_FLT ctol, const SX_INT verb)
 *
 * \brief Iterative on the coarset level
 *
 * \pars  A         pointer to matrix data
 * \pars  b         pointer to rhs data
 * \pars  x         pointer to sol data
 * \pars  ctol      tolerance for the coarsest level
 * \pars  verb   level of output
 *
 */
void sx_amg_coarest_solve(SX_MAT *A, SX_VEC *b, SX_VEC *x, const SX_FLT ctol, const SX_INT verb)
{
    const SX_INT n = A->num_rows;
    const SX_INT maxit = SX_MAX(250, SX_MIN(n * n, 1000));

    SX_INT status;
    SX_KRYLOV ks;

    /* try cg first */
    ks.A = A;
    ks.b = b;
    ks.u = x;
    ks.tol = ctol;
    ks.maxit = maxit;
    ks.stop_type = 1;
    ks.verb = 0;
    status = sx_solver_cg(&ks);

    /* try GMRES if cg fails */
    if (status < 0) {
        ks.restart = MAX_RESTART;
        status = sx_solver_gmres(&ks);
    }

    if (status < 0 && verb >= 2) {
        sx_printf("### WARNING: Coarse level solver failed to converge!\n");
    }
}
