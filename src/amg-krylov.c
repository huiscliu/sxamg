
#include "amg-krylov.h"
#include "op-blas.h"
#include "internal.h"

extern SX_RTN sx_solver_amg_solve(SX_AMG *mg, SX_VEC *x, SX_VEC *b);

SX_RTN sx_solver_cg_itnl(SX_KRYLOV *ks, SX_AMG *mg)
{
    SX_MAT *A = ks->A;
    SX_VEC *bg = ks->b;
    SX_VEC *xg = ks->u;
    SX_RTN rtn;

    int itr_max = ks->maxit;

    double tol_rel = ks->tol;
    double b_norm;

    SX_VEC zg;
    SX_VEC rg;
    SX_VEC qg;
    SX_VEC pg;

    double rho0 = 0, rho1 = 0, beta;
    double residual, ires;
    double tol = -1, err_rel = 0;
    int itr_out;
    int n;

    assert(itr_max > 0);
    assert(tol_rel >= 0);

    n = A->num_rows;
    zg = sx_vec_create(n);
    rg = sx_vec_create(n);
    pg = sx_vec_create(n);
    qg = sx_vec_create(n);

    b_norm = sx_blas_vec_norm2(bg);

    /* init */
    sx_blas_mv_amxpbyz(-1, A, xg, 1, bg, &rg);
    ires = residual = sx_blas_vec_norm2(&rg);

    if (residual == 0.) {
        itr_out = 0;
        rtn.ares = 0;
        rtn.rres = 0;
        rtn.nits = 0;

        goto end;
    }

    err_rel = residual;
    tol = tol_rel * err_rel;

    for (itr_out = 0; itr_out < itr_max; itr_out++) {
        double alpha;

        if (mg == NULL) {
            sx_vec_cp(&rg, &zg);
        }
        else {
            sx_blas_vec_set(&zg, 0);
            sx_solver_amg_solve(mg, &zg, &rg);
        }

        rho1 = sx_blas_vec_dot(&zg, &rg);

        if (itr_out == 0) {
            sx_vec_cp(&zg, &pg);
        }
        else {
            beta = rho1 / rho0;

            sx_blas_vec_axpby(1, &zg, beta, &pg);
        }

        sx_blas_mv_mxy(A, &pg, &qg);
        alpha = sx_blas_vec_dot(&qg, &pg);

        alpha = rho1 / alpha;
        rho0 = rho1;

        sx_blas_vec_axpby(alpha, &pg, 1, xg);
        sx_blas_vec_axpby(-alpha, &qg, 1, &rg);

        residual = sx_blas_vec_norm2(&rg);

        if (ks->verb >= 1) {
            printf("cg: itr: %5d, abs res: %.6e, rel res: %.6e, rbn: %.6e\n",
                 itr_out, residual, (err_rel == 0 ? 0 : residual / err_rel),
                 (b_norm == 0 ? 0 : residual / b_norm));
        }

        if (residual <= tol) break;
    }

    if (itr_out < itr_max) itr_out += 1;

    /* out */
    rtn.ares = residual;
    rtn.rres = residual / ires;
    rtn.nits = itr_out;

 end:

    /* out */
    ks->rtn = rtn;

    sx_vec_destroy(&zg);
    sx_vec_destroy(&rg);
    sx_vec_destroy(&pg);
    sx_vec_destroy(&qg);

    return rtn;
}

SX_RTN sx_solver_gmres_itnl(SX_KRYLOV *ks, SX_AMG *mg)
{
    SX_MAT *A = ks->A;
    SX_VEC *rhs = ks->b;
    SX_VEC *x = ks->u;
    SX_FLT tol_rel = ks->tol;
    SX_INT itr_max = ks->maxit;
    SX_INT m = ks->restart;
    SX_INT verb = ks->verb;
    SX_RTN rtn;

    double b_norm;

    SX_VEC wj;
    SX_VEC *v;

    double **Hg;
    double *ym;
    double *gg;
    SX_VEC rg;
    double *c, *s;

    int itr_inner;
    int i, j = 0;
    int n;
    double tol = 0, err_rel = 0;
    double beta, ires;

    assert(itr_max >= 0);
    assert(ks->tol >= 0);

    n = A->num_rows;
    wj = sx_vec_create(n);
    rg = sx_vec_create(n);

    v = malloc(m * sizeof(*v));
    for (i = 0; i < m; i++) {
        v[i] = sx_vec_create(n);
    }

    gg = malloc((m + 1) * sizeof(*gg));
    ym = malloc(m * sizeof(*ym));

    Hg = (double **)malloc((m + 1) * sizeof(double *));
    for (i = 0; i <= m; i++) {
        Hg[i] = (double *)malloc(m * sizeof(double));
    }

    c = (double *) malloc(m * sizeof(*c));
    s = (double *) malloc(m * sizeof(*s));

    b_norm = sx_blas_vec_norm2(rhs);

    /* init */
    sx_blas_mv_amxpbyz(-1, A, x, 1, rhs, &rg);
    ires = beta = sx_blas_vec_norm2(&rg);

    if (beta == 0) {
        itr_inner = 0;
        rtn.ares = rtn.rres = 0;
        rtn.nits = 0;
        goto end;
    }

    err_rel = beta;
    tol = tol_rel * err_rel;

    itr_inner = 0;
    while (itr_inner < itr_max) {
        int kk;
        double h1, h2, gma;

        /* iit */
        for (kk = 1; kk <= m; kk++) {
            gg[kk] = 0;
        }

        for (kk = 0; kk <= m; kk++) {
            for (i = 0; i < m; i++) Hg[kk][i] = 0;
        }

        /* ortho */
        gg[0] = beta = sx_blas_vec_norm2(&rg);
        sx_blas_vec_axy(1 / beta, &rg, &v[0]);

        for (i = 0; i < m && itr_inner < itr_max; i++) {
            double hij;

            itr_inner++;

            if (mg == NULL) {
                sx_vec_cp(&v[i], &rg);
            }
            else {
                sx_blas_vec_set(&rg, 0);
                sx_solver_amg_solve(mg, &rg, &v[i]);
            }

            sx_blas_mv_mxy(A, &rg, &wj);

            for (j = 0; j <= i; j++) {
                hij = sx_blas_vec_dot(&wj, &v[j]);
                sx_blas_vec_axpby(-hij, &v[j], 1, &wj);

                Hg[j][i] = hij;
            }

            hij = sx_blas_vec_norm2(&wj);
            Hg[i + 1][i] = hij;

            /* break down */
            if (fabs(hij) <= 1e-30) {
                i -= 1;
                break;
            }
            else if (i + 1 < m) {
                sx_blas_vec_axy(1 / hij, &wj, &v[i + 1]);
            }

            for (j = 0; j < i; j++) {
                h1 = c[j] * Hg[j][i] + s[j] * Hg[j + 1][i];
                h2 = -s[j] * Hg[j][i] + c[j] * Hg[j + 1][i];

                Hg[j][i] = h1;
                Hg[j + 1][i] = h2;
            }

            gma = sqrt(Hg[i][i] * Hg[i][i] + Hg[i + 1][i] * Hg[i + 1][i]);
            if (fabs(gma) == 0.) gma = 1e-20;

            c[i] = Hg[i][i] / gma;
            s[i] = Hg[i + 1][i] / gma;

            gg[i + 1] = -s[i] * gg[i];
            gg[i] = c[i] * gg[i];
            Hg[i][i] = c[i] * Hg[i][i] + s[i] * Hg[i + 1][i];
            beta = fabs(gg[i + 1]);

            if (verb >= 1) {
                printf("gmres: itr: %4d, abs res: %.6e, rel res: %.6e, rbn: %.6e\n", \
                        itr_inner, beta, (err_rel == 0 ? 0 : beta / err_rel), \
                        (b_norm == 0 ? 0 : beta / b_norm));
            }

            if (beta <= tol) {
                goto solve;
            }
        }

solve:
        kk = i == m ? m : i + 1;
        for (i = kk - 1; i >= 0; i--) {
            ym[i] = gg[i] / Hg[i][i];
            for (j = 0; j < i; j++) {
                gg[j] = gg[j] - ym[i] * Hg[j][i];
            }
        }

        /* Update z */
        if (kk > 0) {
            sx_blas_vec_axy(ym[kk - 1], &v[kk - 1], &rg);

            for (i = kk - 2; i >= 0; i--) {
                sx_blas_vec_axpby(ym[i], &v[i], 1, &rg);
            }

            if (mg == NULL) {
                sx_blas_vec_axpby(1, &rg, 1, x);
            }
            else {
                sx_blas_vec_set(&wj, 0);
                sx_solver_amg_solve(mg, &wj, &rg);
                sx_blas_vec_axpby(1, &wj, 1, x);
            }
        }

        if (beta <= tol) {
            break;
        }
        else {
            sx_blas_mv_amxpbyz(-1, A, x, 1, rhs, &rg);
        }
    }

    /* out */
    rtn.ares = beta;
    rtn.rres = beta / ires;
    rtn.nits = itr_inner;

end:

    /* out */
    ks->rtn = rtn;

    sx_vec_destroy(&wj);
    sx_vec_destroy(&rg);

    for (i = 0; i < m; i++) {
        sx_vec_destroy(&v[i]);
    }
    free(v);

    free(c);
    free(s);
    free(gg);
    free(ym);

    for (i = 0; i <= m; i++) {
        free(Hg[i]);
    }
    free(Hg);

    return rtn;
}

SX_RTN sx_solver_bicgstab_itnl(SX_KRYLOV *ks, SX_AMG *mg)
{
    SX_MAT *A = ks->A;
    SX_VEC *bg = ks->b;
    SX_VEC *xg = ks->u;
    SX_RTN rtn;

    int itr_max = ks->maxit;
    double tol_rel = ks->tol;
    double b_norm;

    SX_VEC rg;
    SX_VEC rh;
    SX_VEC pg;
    SX_VEC ph;
    SX_VEC sg;
    SX_VEC sh;
    SX_VEC tg;
    SX_VEC vg;

    double rho0 = 0, rho1 = 0;
    double alpha = 0, beta = 0, omega = 0;
    double residual;
    double tol = -1, err_rel = 0;
    int itr_out;
    int i;
    int n;

    assert(itr_max >= 0);

    n = A->num_rows;

    rg = sx_vec_create(n);
    rh = sx_vec_create(n);
    pg = sx_vec_create(n);
    ph = sx_vec_create(n);
    sg = sx_vec_create(n);
    sh = sx_vec_create(n);
    tg = sx_vec_create(n);
    vg = sx_vec_create(n);

    /* init */
    sx_blas_mv_amxpbyz(-1, A, xg, 1, bg, &rg);
    sx_vec_cp(&rg, &rh);
    sx_blas_vec_set(&sh, 0);
    sx_blas_vec_set(&ph, 0);

    b_norm = sx_blas_vec_norm2(bg);

    residual = err_rel = sx_blas_vec_norm2(&rg);
    if (residual == 0.) {
        itr_out = 0;
        rtn.ares = rtn.rres = 0;
        rtn.nits = 0;
        goto end;
    }

    tol = residual * tol_rel;

    for (itr_out = 0; itr_out < itr_max; itr_out++) {
        rho1 = sx_blas_vec_dot(&rg, &rh);

        if (rho1 == 0) {
            printf("bicgstab: method failed.!\n");
            break;
        }

        if (itr_out == 0) {
            for (i = 0; i < n; i++) pg.d[i] = rg.d[i];
        }
        else {
            beta = (rho1 * alpha) / (rho0 * omega);
            for (i = 0; i < n; i++) {
                pg.d[i] = rg.d[i] + beta * (pg.d[i] - omega * vg.d[i]);
            }
        }

        rho0 = rho1;

        if (mg == NULL) {
            sx_vec_cp(&pg, &ph);
        }
        else {
            sx_vec_set_value(&ph, 0.);
            sx_solver_amg_solve(mg, &ph, &pg);
        }

        sx_blas_mv_amxpbyz(1, A, &ph, 0, &pg, &vg);

        alpha = rho1 / sx_blas_vec_dot(&rh, &vg);

        sx_blas_vec_axpbyz(1, &rg, -alpha, &vg, &sg);

        if (sx_blas_vec_norm2(&sg) <= 1e-30) {
            printf("bicgstab: ||s|| is too small: %f, terminated.\n", sx_blas_vec_norm2(&sg));

            sx_blas_vec_axpby(alpha, &ph, 1, xg);

            sx_blas_mv_amxpbyz(-1, A, xg, 1, bg, &rg);
            residual = sx_blas_vec_norm2(&rg);

            break;
        }

        if (mg == NULL) {
            sx_vec_cp(&sg, &sh);
        }
        else {
            sx_vec_set_value(&sh, 0.);
            sx_solver_amg_solve(mg, &sh, &sg);
        }

        sx_blas_mv_amxpbyz(1, A, &sh, 0, &pg, &tg);

        omega = sx_blas_vec_dot(&tg, &sg) / sx_blas_vec_dot(&tg, &tg);
        for (i = 0; i < n; i++) {
            xg->d[i] = xg->d[i] + alpha * ph.d[i] + omega * sh.d[i];
            rg.d[i] = sg.d[i] - omega * tg.d[i];
        }

        residual = sx_blas_vec_norm2(&rg);

        if (ks->verb >= 1) {
            printf("bicgstab: itr: %5d, abs res: %.6e, rel res: %.6e, "
                    "rbn: %.6e\n", itr_out, residual, (err_rel == 0 ? 0 : residual / err_rel), 
                    (b_norm == 0 ? 0 : residual / b_norm));
        }

        if (residual <= tol) break;
    }

    if (itr_out < itr_max) itr_out += 1;

    rtn.ares = residual;
    rtn.rres = residual / err_rel;
    rtn.nits = itr_out;

end:

    ks->rtn = rtn;

    sx_vec_destroy(&rg);
    sx_vec_destroy(&rh);
    sx_vec_destroy(&pg);
    sx_vec_destroy(&ph);
    sx_vec_destroy(&sg);
    sx_vec_destroy(&sh);
    sx_vec_destroy(&tg);
    sx_vec_destroy(&vg);

    return rtn;
}
