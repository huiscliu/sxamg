
#include "sxamg.h"
#include "mat.h"

int main(void)
{
    SX_MAT A;
    SX_INT ncx = 23, ncy = 13, ncz = 24;
    SX_INT prob = 3;
    SX_INT nglobal = 0;
    SX_INT k, i;
    SX_VEC x, b;

    assert(ncx > 0);
    assert(ncy > 0);
    assert(ncz > 0);

    /* create distributed matrix */
    if (prob == 1) {
        nglobal = ncx;
        A = laplacian_3pt(ncx);
    }
    else if (prob == 2) {
        nglobal = ncx * ncx;
        A = laplacian_5pt(ncx);
    }
    else if (prob == 3) {
        nglobal = ncx * ncy * ncz;
        A = laplacian_7pt_bound(ncx, ncy, ncz);
    }
    else {
        sx_printf("sx: wrong problem.\n");
        exit(-1);
    }

    /* create vectors */
    k = nglobal;
    x = sx_vec_create(nglobal);
    b = sx_vec_create(nglobal);

    /* init x and b */
    for (i = 0; i < k; i++) {
        /* right hand size */
        sx_vec_set_entry(&b, i, 1);

        /* initial guess */
        sx_vec_set_entry(&x, i, 0);
    }

    /* solve Ax = b using CG with AMG preconditioner */
    {
        SX_INT maxits = 2000;  /* maximal iterations */
        SX_FLT tol = 1e-4;     /* stop tolerance */
        SX_FLT err, err0;
        SX_VEC r;     /* residual */
        SX_VEC p;
        SX_VEC z, q;
        SX_FLT rho1, rho0 = 0., alpha, beta;

        /* preconditioner */
        SX_AMG mg;
        SX_AMG_PARS pars;

        /* pars */
        sx_amg_pars_init(&pars);
        pars.maxit = 1;
        sx_amg_setup(&mg, &A, &pars);

        /* create vector */
        r = sx_vec_create(k);
        z = sx_vec_create(k);
        q = sx_vec_create(k);
        p = sx_vec_create(k);

        /* initial residual */
        sx_blas_mv_amxpbyz(-1., &A, &x, 1., &b, &r);
        err0 = sx_blas_vec_norm2(&r);

        sx_printf("Convergence settings: relative residual: %"fFMTe
                ", maximal iterations: %"dFMT"\n\n", tol, maxits);

        sx_printf("Initial residual: %"fFMTe"\n\n", err0);

        for (i = 0; i < maxits; i++) {
            /* supposed to solve preconditioning system Mz = r */
#if 0       /* no pc */
            sx_blas_vec_copy(&r, &z);
#else       /* amg pc */
            sx_blas_vec_set(&z, 0.);
#endif
            sx_solver_amg_solve(&mg, &z, &r);

            /* rho = <r, z> */
            rho1 = sx_blas_vec_dot(&r, &z);

            if (i == 0) {
                /* p = z */
                sx_blas_vec_copy(&z, &p);
            }
            else {
                beta = rho1 / rho0;

                /* update p */
                sx_blas_vec_axpby(1, &z, beta, &p);
            }

            /* save rho */
            rho0 = rho1;

            /* update q */
            sx_blas_mv_mxy(&A, &p, &q);

            /* compute alpha */
            alpha = rho1 / sx_blas_vec_dot(&p, &q);

            /* update x */
            sx_blas_vec_axpy(alpha, &p, &x);

            /* update r */
            sx_blas_vec_axpy(-alpha, &q, &r);

            /* check convergence */
            err = sx_blas_vec_norm2(&r);

            sx_printf("itr: %6"dFMT",     residual: %"fFMTe", relative error: %"fFMTe"\n",
                    i + 1, err, err / err0);

            if (err / err0 <= tol) break;
        }

        if (i == maxits) {
            sx_printf("\ncg failed to converge.\n\n");
        }
        else {
            sx_printf("\ncg converged: absolute residual: %"fFMTe
                    ", total iterations: %"dFMT"\n\n", err, i + 1);
        }

        sx_vec_destroy(&r);
        sx_vec_destroy(&p);
        sx_vec_destroy(&z);
        sx_vec_destroy(&q);

        sx_amg_data_destroy(&mg);
    }

    /* release memory */
    sx_mat_destroy(&A);
    sx_vec_destroy(&x);
    sx_vec_destroy(&b);

    return 0;
}
