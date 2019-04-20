
#include "op-smoother.h"
#include "internal.h"
#include <assert.h>

static void Diaginv(SX_MAT *, SX_FLT *);
static SX_FLT DinvAnorminf(SX_MAT *, SX_FLT *);
static void Diagx(SX_FLT *, SX_INT, SX_FLT *, SX_FLT *);

static void Rr(SX_MAT *, SX_FLT *, SX_FLT *, SX_FLT *, SX_FLT *, SX_FLT *,
        SX_FLT *, SX_FLT *, SX_INT);

/**
 * \fn void sx_amg_smoother_poly (SX_MAT *Amat, SX_VEC *brhs, SX_VEC *usol, 
 *                                   SX_INT n, SX_INT ndeg, SX_INT L)
 *
 * \brief poly approx to A^{-1} as MG smoother 
 *
 * \pars Amat  Pointer to stiffness matrix, consider square matrix.
 * \pars brhs  Pointer to right hand side
 * \pars usol  Pointer to solution 
 * \pars n     Problem size 
 * \pars ndeg  Degree of poly 
 * \pars L     Number of iterations
 *
 */
static void sx_amg_smoother_poly(SX_MAT * Amat, SX_VEC * brhs, SX_VEC * usol, SX_INT n,
        SX_INT ndeg, SX_INT L)
{
    // local variables
    SX_INT i;
    SX_FLT *b = brhs->d, *u = usol->d;
    SX_FLT *Dinv = NULL, *r = NULL, *rbar = NULL, *v0 = NULL, *v1 = NULL;
    SX_FLT *error = NULL, *k = NULL;
    SX_FLT mu0, mu1, smu0, smu1;

    /* allocate memory */
    Dinv = (SX_FLT *) sx_calloc(n, sizeof(SX_FLT));
    r = (SX_FLT *) sx_calloc(n, sizeof(SX_FLT));
    rbar = (SX_FLT *) sx_calloc(n, sizeof(SX_FLT));
    v0 = (SX_FLT *) sx_calloc(n, sizeof(SX_FLT));
    v1 = (SX_FLT *) sx_calloc(n, sizeof(SX_FLT));
    error = (SX_FLT *) sx_calloc(n, sizeof(SX_FLT));
    k = (SX_FLT *) sx_calloc(6, sizeof(SX_FLT));    // coefficients for calculation

    // get the inverse of the diagonal of A
    Diaginv(Amat, Dinv);

    // set up parseter
    mu0 = DinvAnorminf(Amat, Dinv);     // get the inf norm of Dinv*A;

    mu0 = 1.0 / mu0;
    mu1 = 4.0 * mu0;            // default set 8;
    smu0 = sqrt(mu0);
    smu1 = sqrt(mu1);

    k[1] = (mu0 + mu1) / 2.0;
    k[2] = (smu0 + smu1) * (smu0 + smu1) / 2.0;
    k[3] = mu0 * mu1;

    // 4.0*mu0*mu1/(sqrt(mu0)+sqrt(mu1))/(sqrt(mu0)+sqrt(mu1));
    k[4] = 2.0 * k[3] / k[2];

    // square of (sqrt(kappa)-1)/(sqrt(kappa)+1);
    k[5] = (mu1 - 2.0 * smu0 * smu1 + mu0) / (mu1 + 2.0 * smu0 * smu1 + mu0);

    // Update 
    for (i = 0; i < L; i++) {
        SX_VEC vr;

        vr.n = usol->n;
        vr.d = r;

        // get residual
        sx_blas_mv_mxy(Amat, usol, &vr); // r= Amat*u;
        sx_blas_array_axpyz(n, -1, r, b, r);  // r= -r+b;

        // Get correction error = R*r
        Rr(Amat, Dinv, r, rbar, v0, v1, error, k, ndeg);

        // update solution
        sx_blas_array_axpy(n, 1, error, u);

    }

    // free memory
    sx_free(Dinv);
    sx_free(r);
    sx_free(rbar);
    sx_free(v0);
    sx_free(v1);
    sx_free(error);
    sx_free(k);

    return;

}

/**
 * \fn static void Diaginv (SX_MAT *Amat, SX_FLT *Dinv)
 *
 * \brief find the inverse of the diagonal of A
 * 
 * \pars Amat sparse matrix in CSR format
 * \pars Dinv vector to store the inverse of diagonal of A
 *
 */
static void Diaginv(SX_MAT * Amat, SX_FLT * Dinv)
{
    const SX_INT n = Amat->num_rows;
    const SX_INT *ia = Amat->Ap, *ja = Amat->Aj;
    const SX_FLT *a = Amat->Ax;
    SX_INT i, j;

    for (i = 0; i < n; i++) {
        for (j = ia[i]; j < ia[i + 1]; j++) {
            if (i == ja[j])     // find the diagonal 
                break;
        }

        Dinv[i] = 1.0 / a[j];
    }
}

/**
 * \fn    static SX_FLT DinvAnorminf(SX_MAT *Amat, SX_FLT *Dinv)
 *
 * \brief Get the inf norm of Dinv*A
 * 
 * \pars Amat sparse matrix in CSR format
 * \pars Dinv vector to store the inverse of diagonal of A
 *
 * \return Inf norm of Dinv*A
 *
 */
static SX_FLT DinvAnorminf(SX_MAT * Amat, SX_FLT * Dinv)
{
    //local variable
    const SX_INT n = Amat->num_rows;
    const SX_INT *ia = Amat->Ap;
    const SX_FLT *a = Amat->Ax;
    int i, j;
    SX_FLT norm, temp;

    // get the infinity norm of Dinv*A
    norm = 0.;
    for (i = 0; i < n; i++) {
        temp = 0.;
        for (j = ia[i]; j < ia[i + 1]; j++) {
            temp += SX_ABS(a[j]);
        }
        temp *= Dinv[i];        // temp is the L1 norm of the ith row of Dinv*A;
        norm = SX_MAX(norm, temp);
    }

    return norm;
}

/**
 * \fn    static void Diagx(SX_FLT *Dinv, SX_INT n, SX_FLT *x, SX_FLT *b)
 *
 * \brief Compute b = Dinv * x;
 * 
 * \pars Dinv  Vector to represent Diagonal matrix
 * \pars n     Problem size
 * \pars b     Vector
 * \pars x     Vector 
 *
 */
static void Diagx(SX_FLT * Dinv, SX_INT n, SX_FLT * x, SX_FLT * b)
{
    SX_INT i;

    for (i = 0; i < n; i++) {
        b[i] = Dinv[i] * x[i];
    }
}

/**
 * \fn static void Rr (SX_MAT *Amat, SX_FLT *Dinv, SX_FLT *r, SX_FLT *rbar, 
 *                     SX_FLT *v0, SX_FLT *v1, SX_FLT *vnew, SX_FLT *k, SX_INT m)
 *
 * \brief Compute action R*r, where R =  q_m(Dinv*A)*Dinv;
 * 
 * \pars Amat sparse matrix A in CSR format
 * \pars Dinv inverse of the diagoal of A
 * \pars r  initial residual
 * \pars rbar intermediate memory
 * \pars v0 0 degree action R*r
 * \pars v1 1 degree action R*r
 * \pars vnew store final action vnew = R*r;
 * \pars k  coefficients in equations
 * \pars m  degree of polynomial
 *
 */
static void Rr(SX_MAT * Amat, SX_FLT * Dinv, SX_FLT * r, SX_FLT * rbar, SX_FLT * v0,
        SX_FLT * v1, SX_FLT *vnew, SX_FLT * k, SX_INT m)
{
    // local variables
    const SX_INT n = Amat->num_rows;
    SX_INT i, j;
    SX_VEC vr, vv;

    vr.n = vv.n = n;
    vr.d = rbar;
    vv.d = v1;

    //1 set up rbar
    Diagx(Dinv, n, r, rbar);    // rbar = Dinv *r;

    //2 set up v0, v1;
    sx_blas_mv_mxy(Amat, &vr, &vv); //v1= A*rbar;
    Diagx(Dinv, n, v1, v1);     // v1=Dinv *v1;

    for (i = 0; i < n; i++) {
        v0[i] = k[1] * rbar[i];
        v1[i] = k[2] * rbar[i] - k[3] * v1[i];
    }

    //3 iterate to get v_(j+1)
    for (j = 1; j < m; j++) {
        sx_blas_mv_mxy(Amat, &vv, &vr);     //rbar= A*v_(j);

        for (i = 0; i < n; i++) {
            rbar[i] = (r[i] - rbar[i]) * Dinv[i];       // indeed rbar=Dinv*(r-A*v_(j));
            vnew[i] = v1[i] + k[5] * (v1[i] - v0[i]) + k[4] * rbar[i];  // compute v_(j+1)
            // prepare for next cycle
            v0[i] = v1[i];
            v1[i] = vnew[i];
        }
    }
}

/**
 * \fn void sx_amg_smoother_jacobi (SX_VEC *u, const SX_INT i_1, const SX_INT i_n,
 *                                     const SX_INT s, SX_MAT *A, SX_VEC *b, SX_INT L)
 *
 * \brief Jacobi method as a smoother
 *
 * \pars u      Pointer to SX_VEC: the unknowns (IN: initial, OUT: approximation)
 * \pars i_1    Starting index
 * \pars i_n    Ending index
 * \pars s      Increasing step
 * \pars A      Pointer to dBSRmat: the coefficient matrix
 * \pars b      Pointer to SX_VEC: the right hand side
 * \pars L      Number of iterations
 *
 */
static void sx_amg_smoother_jacobi(SX_VEC * u, const SX_INT i_1, const SX_INT i_n, const SX_INT s,
        SX_MAT * A, SX_VEC * b, SX_INT L)
{
    const SX_INT N = SX_ABS(i_n - i_1) + 1;
    const SX_INT *ia = A->Ap, *ja = A->Aj;
    const SX_FLT *aj = A->Ax, *bval = b->d;
    SX_FLT *uval = u->d;
    SX_FLT w = 0.8;              //0.8
    // local variables
    SX_INT i, j, k, begin_row, end_row;

    SX_FLT *t = (SX_FLT *) sx_calloc(N, sizeof(SX_FLT));
    SX_FLT *d = (SX_FLT *) sx_calloc(N, sizeof(SX_FLT));

    while (L--) {
        if (s > 0) {
            for (i = i_1; i <= i_n; i += s) {
                t[i] = bval[i];
                begin_row = ia[i], end_row = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) {
                    j = ja[k];
                    if (i != j)
                        t[i] -= aj[k] * uval[j];
                    else
                        d[i] = aj[k];
                }
            }
            for (i = i_1; i <= i_n; i += s) {
                if (SX_ABS(d[i]) > SMALLFLOAT)
                    uval[i] = (1 - w) * uval[i] + w * t[i] / d[i];
            }
        }
        else {
            for (i = i_1; i >= i_n; i += s) {
                t[i] = bval[i];
                begin_row = ia[i], end_row = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) {
                    j = ja[k];
                    if (i != j)
                        t[i] -= aj[k] * uval[j];
                    else
                        d[i] = aj[k];
                }
            }

            for (i = i_1; i >= i_n; i += s) {
                //if (ABS(d[i])>SMALLFLOAT) uval[i]=t[i]/d[i];
                if (SX_ABS(d[i]) > SMALLFLOAT)
                    uval[i] = (1 - w) * uval[i] + w * t[i] / d[i];
            }
        }
    }

    sx_free(t);
    sx_free(d);
}

/**
 * \fn void sx_amg_smoother_gs (SX_VEC *u, const SX_INT i_1, const SX_INT i_n,
 *                                 const SX_INT s, SX_MAT *A, SX_VEC *b, SX_INT L)
 *
 * \brief Gauss-Seidel method as a smoother
 *
 * \pars u    Pointer to SX_VEC: the unknowns (IN: initial, OUT: approximation)
 * \pars i_1  Starting index
 * \pars i_n  Ending index
 * \pars s    Increasing step
 * \pars A    Pointer to dBSRmat: the coefficient matrix
 * \pars b    Pointer to SX_VEC: the right hand side
 * \pars L    Number of iterations
 *
 */
static void sx_amg_smoother_gs(SX_VEC * u, const SX_INT i_1, const SX_INT i_n, const SX_INT s,
        SX_MAT *A, SX_VEC *b, SX_INT L)
{
    const SX_INT *ia = A->Ap, *ja = A->Aj;
    const SX_FLT *aj = A->Ax, *bval = b->d;
    SX_FLT *uval = u->d;

    // local variables
    SX_INT i, j, k, begin_row, end_row;
    SX_FLT t, d = 0.0;

    if (s > 0) {
        while (L--) {
            for (i = i_1; i <= i_n; i += s) {
                t = bval[i];
                begin_row = ia[i], end_row = ia[i + 1];

                for (k = begin_row; k < end_row; ++k) {
                    j = ja[k];
                    if (i != j)
                        t -= aj[k] * uval[j];
                    else if (SX_ABS(aj[k]) > SMALLFLOAT)
                        d = 1.e+0 / aj[k];
                }

                uval[i] = t * d;
            }                   // end for i
        }                       // end while
    }                           // if s
    else {

        while (L--) {
            for (i = i_1; i >= i_n; i += s) {
                t = bval[i];
                begin_row = ia[i], end_row = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) {
                    j = ja[k];
                    if (i != j)
                        t -= aj[k] * uval[j];
                    else if (SX_ABS(aj[k]) > SMALLFLOAT)
                        d = 1.0 / aj[k];
                }

                uval[i] = t * d;
            }                   // end for i
        }                       // end while
    }                           // end if
}

/**
 * \fn void sx_amg_smoother_gs_cf (SX_VEC *u, SX_MAT *A, SX_VEC *b, SX_INT L,
 *                                    SX_INT *mark, const SX_INT order)
 *
 * \brief Gauss-Seidel smoother with C/F ordering for Au=b
 *
 * \pars u      Pointer to SX_VEC: the unknowns (IN: initial, OUT: approximation)
 * \pars A      Pointer to dBSRmat: the coefficient matrix
 * \pars b      Pointer to SX_VEC: the right hand side
 * \pars L      Number of iterations
 * \pars mark   C/F marker array
 * \pars order  C/F ordering: -1: F-first; 1: C-first
 *
 */
static void sx_amg_smoother_gs_cf(SX_VEC * u, SX_MAT * A, SX_VEC * b, SX_INT L, SX_INT * mark,
        const SX_INT order)
{
    const SX_INT nrow = b->n;    // number of rows
    const SX_INT *ia = A->Ap, *ja = A->Aj;
    const SX_FLT *aj = A->Ax, *bval = b->d;
    SX_FLT *uval = u->d;

    SX_INT i, j, k, begin_row, end_row;
    SX_FLT t, d = 0.0;

    // F-point first, C-point second
    if (order) {
        while (L--) {
            for (i = 0; i < nrow; i++) {
                if (mark[i] != 1) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i + 1];
                    for (k = begin_row; k < end_row; k++) {
                        j = ja[k];
                        if (i != j)
                            t -= aj[k] * uval[j];
                        else
                            d = aj[k];
                    }           // end for k

                    if (SX_ABS(d) > SMALLFLOAT) uval[i] = t / d;
                }
            }

            for (i = 0; i < nrow; i++) {
                if (mark[i] == 1) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i + 1];
                    for (k = begin_row; k < end_row; k++) {
                        j = ja[k];
                        if (i != j)
                            t -= aj[k] * uval[j];
                        else
                            d = aj[k];
                    }           // end for k

                    if (SX_ABS(d) > SMALLFLOAT) uval[i] = t / d;
                }
            }                   // end for i
        }                       // end while
    }
    else {                      // C-point first, F-point second
        while (L--) {
            for (i = 0; i < nrow; i++) {
                if (mark[i] == 1) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i + 1];
                    for (k = begin_row; k < end_row; k++) {
                        j = ja[k];
                        if (i != j)
                            t -= aj[k] * uval[j];
                        else
                            d = aj[k];
                    }           // end for k

                    if (SX_ABS(d) > SMALLFLOAT) uval[i] = t / d;
                }
            }

            for (i = 0; i < nrow; i++) {
                if (mark[i] != 1) {
                    t = bval[i];
                    begin_row = ia[i], end_row = ia[i + 1];
                    for (k = begin_row; k < end_row; k++) {
                        j = ja[k];
                        if (i != j)
                            t -= aj[k] * uval[j];
                        else
                            d = aj[k];
                    }           // end for k

                    if (SX_ABS(d) > SMALLFLOAT)
                        uval[i] = t / d;
                }
            }                   // end for i
        }                       // end while
    }                           // end if order
}

/**
 * \fn void sx_amg_smoother_sgs (SX_VEC *u, SX_MAT *A, SX_VEC *b, SX_INT L)
 *
 * \brief Symmetric Gauss-Seidel method as a smoother
 *
 * \pars u      Pointer to SX_VEC: the unknowns (IN: initial, OUT: approximation)
 * \pars A      Pointer to dBSRmat: the coefficient matrix
 * \pars b      Pointer to SX_VEC: the right hand side
 * \pars L      Number of iterations
 *
 */
static void sx_amg_smoother_sgs(SX_VEC * u, SX_MAT * A, SX_VEC * b, SX_INT L)
{
    const SX_INT nm1 = b->n - 1;
    const SX_INT *ia = A->Ap, *ja = A->Aj;
    const SX_FLT *aj = A->Ax, *bval = b->d;
    SX_FLT *uval = u->d;

    // local variables
    SX_INT i, j, k, begin_row, end_row;
    SX_FLT t, d = 0;

    while (L--) {
        // forward sweep
        for (i = 0; i <= nm1; ++i) {
            t = bval[i];
            begin_row = ia[i], end_row = ia[i + 1];
            for (k = begin_row; k < end_row; ++k) {
                j = ja[k];
                if (i != j)
                    t -= aj[k] * uval[j];
                else
                    d = aj[k];
            }                   // end for k

            if (SX_ABS(d) > SMALLFLOAT) uval[i] = t / d;
        }                       // end for i

        // backward sweep
        for (i = nm1 - 1; i >= 0; --i) {
            t = bval[i];
            begin_row = ia[i], end_row = ia[i + 1];
            for (k = begin_row; k < end_row; ++k) {
                j = ja[k];
                if (i != j)
                    t -= aj[k] * uval[j];
                else
                    d = aj[k];
            }                   // end for k

            if (SX_ABS(d) > SMALLFLOAT) uval[i] = t / d;
        }
    }
}

/**
 * \fn void sx_amg_smoother_sor (SX_VEC *u, const SX_INT i_1, const SX_INT i_n, const SX_INT s,
 *                                  SX_MAT *A, SX_VEC *b, SX_INT L, const SX_FLT w)
 *
 * \brief SOR method as a smoother
 *
 * \pars u      Pointer to SX_VEC: the unknowns (IN: initial, OUT: approximation)
 * \pars i_1    Starting index
 * \pars i_n    Ending index
 * \pars s      Increasing step
 * \pars A      Pointer to dBSRmat: the coefficient matrix
 * \pars b      Pointer to SX_VEC: the right hand side
 * \pars L      Number of iterations
 * \pars w      Over-relaxation weight
 *
 */
static void sx_amg_smoother_sor(SX_VEC * u, const SX_INT i_1, const SX_INT i_n, const SX_INT s,
        SX_MAT *A, SX_VEC *b, SX_INT L, const SX_FLT w)
{
    const SX_INT *ia = A->Ap, *ja = A->Aj;
    const SX_FLT *aj = A->Ax, *bval = b->d;
    SX_FLT *uval = u->d;

    // local variables
    SX_INT i, j, k, begin_row, end_row;
    SX_FLT t, d = 0;

    while (L--) {
        if (s > 0) {
            for (i = i_1; i <= i_n; i += s) {
                t = bval[i];
                begin_row = ia[i], end_row = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) {
                    j = ja[k];
                    if (i != j)
                        t -= aj[k] * uval[j];
                    else
                        d = aj[k];
                }

                if (SX_ABS(d) > SMALLFLOAT)
                    uval[i] = w * (t / d) + (1 - w) * uval[i];
            }
        }
        else {
            for (i = i_1; i >= i_n; i += s) {
                t = bval[i];
                begin_row = ia[i], end_row = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) {
                    j = ja[k];
                    if (i != j)
                        t -= aj[k] * uval[j];
                    else
                        d = aj[k];
                }

                if (SX_ABS(d) > SMALLFLOAT)
                    uval[i] = w * (t / d) + (1 - w) * uval[i];
            }
        }
    }
}

/**
 * \fn void sx_amg_smoother_L1diag (SX_VEC *u, const SX_INT i_1, const SX_INT i_n, const SX_INT s,
 *                                     SX_MAT *A, SX_VEC *b, SX_INT L)
 *
 * \brief Diagonal scaling (using L1 norm) as a smoother
 *
 * \pars u      Pointer to SX_VEC: the unknowns (IN: initial, OUT: approximation)
 * \pars i_1    Starting index
 * \pars i_n    Ending index
 * \pars s      Increasing step
 * \pars A      Pointer to dBSRmat: the coefficient matrix
 * \pars b      Pointer to SX_VEC: the right hand side
 * \pars L      Number of iterations
 *
 */
static void sx_amg_smoother_L1diag(SX_VEC * u, const SX_INT i_1, const SX_INT i_n, const SX_INT s,
        SX_MAT *A, SX_VEC *b, SX_INT L)
{
    const SX_INT N = SX_ABS(i_n - i_1) + 1;
    const SX_INT *ia = A->Ap, *ja = A->Aj;
    const SX_FLT *aj = A->Ax, *bval = b->d;
    SX_FLT *uval = u->d;

    // local variables
    SX_INT i, j, k, begin_row, end_row;

    // Checks should be outside of for; t,d can be allocated before calling
    SX_FLT *t = (SX_FLT *) sx_calloc(N, sizeof(SX_FLT));
    SX_FLT *d = (SX_FLT *) sx_calloc(N, sizeof(SX_FLT));

    while (L--) {
        if (s > 0) {
            for (i = i_1; i <= i_n; i += s) {
                t[i] = bval[i];
                d[i] = 0.0;
                begin_row = ia[i], end_row = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) {
                    j = ja[k];
                    t[i] -= aj[k] * uval[j];
                    d[i] += SX_ABS(aj[k]);
                }
            }

            for (i = i_1; i <= i_n; i += s) {
                if (SX_ABS(d[i]) > SMALLFLOAT)
                    u->d[i] += t[i] / d[i];
            }
        }
        else {
            for (i = i_1; i >= i_n; i += s) {
                t[i] = bval[i];
                d[i] = 0.0;
                begin_row = ia[i], end_row = ia[i + 1];
                for (k = begin_row; k < end_row; ++k) {
                    j = ja[k];
                    t[i] -= aj[k] * uval[j];
                    d[i] += SX_ABS(aj[k]);
                }
            }

            for (i = i_1; i >= i_n; i += s) {
                if (SX_ABS(d[i]) > SMALLFLOAT) u->d[i] += t[i] / d[i];
            }
        }
    }

    sx_free(t);
    sx_free(d);
}

/**
 * \fn static void sx_amg_smoother_pre (SX_SMTR *s)
 *
 * \brief Multigrid presmoothing
 *
 */
void sx_amg_smoother_pre(SX_SMTR *s)
{
    SX_SM_TYPE smoother;
    SX_MAT * A;
    SX_VEC * b;
    SX_VEC * x;
    SX_INT nsweeps;
    SX_INT istart;
    SX_INT iend;
    SX_INT istep;
    SX_FLT relax;
    SX_INT ndeg;
    SX_INT order;
    SX_INT *ordering;

    assert(s != NULL);

    smoother = s->smoother;
    A = s->A;
    b = s->b;
    x = s->x;
    nsweeps = s->nsweeps;
    istart = s->istart;
    iend = s->iend;
    istep = s->istep;
    relax = s->relax;
    ndeg = s->ndeg;
    order = s->cf_order;
    ordering = s->ordering;

    switch (smoother) {
        case SX_SM_GS:
            if (order && ordering != NULL) {
                sx_amg_smoother_gs_cf(x, A, b, nsweeps, ordering, 1);
            }
            else {
                sx_amg_smoother_gs(x, istart, iend, istep, A, b, nsweeps);
            }
            break;

        case SX_SM_SGS:
            sx_amg_smoother_sgs(x, A, b, nsweeps);
            break;

        case SX_SM_JACOBI:
            sx_amg_smoother_jacobi(x, istart, iend, istep, A, b, nsweeps);
            break;

        case SX_SM_L1DIAG:
            sx_amg_smoother_L1diag(x, istart, iend, istep, A, b, nsweeps);
            break;

        case SX_SM_POLY:
            sx_amg_smoother_poly(A, b, x, iend + 1, ndeg, nsweeps);
            break;

        case SX_SM_SOR:
            sx_amg_smoother_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            break;

        case SX_SM_SSOR:
            sx_amg_smoother_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            sx_amg_smoother_sor(x, iend, istart, -istep, A, b, nsweeps, relax);
            break;

        case SX_SM_GSOR:
            sx_amg_smoother_gs(x, istart, iend, istep, A, b, nsweeps);
            sx_amg_smoother_sor(x, iend, istart, -istep, A, b, nsweeps, relax);
            break;

        case SX_SM_SGSOR:
            sx_amg_smoother_gs(x, istart, iend, istep, A, b, nsweeps);
            sx_amg_smoother_gs(x, iend, istart, -istep, A, b, nsweeps);
            sx_amg_smoother_sor(x, istart, iend, istep, A, b, nsweeps, relax);
            sx_amg_smoother_sor(x, iend, istart, -istep, A, b, nsweeps, relax);
            break;

        default:
            sx_printf("### ERROR: Wrong smoother type %"dFMT"!\n", smoother);
            sx_exit_on_errcode(ERROR_INPUT_PAR, __FUNCTION__);
    }
}

/**
 * \fn static void sx_amg_smoother_post (SX_SMTR *s)
 *
 * \brief Multigrid presmoothing
 *
 *
 */
void sx_amg_smoother_post(SX_SMTR *s)
{
    SX_SM_TYPE smoother;
    SX_MAT * A;
    SX_VEC * b;
    SX_VEC * x;
    SX_INT nsweeps;
    SX_INT istart;
    SX_INT iend;
    SX_INT istep;
    SX_FLT relax;
    SX_INT ndeg;
    SX_INT order;
    SX_INT *ordering;

    assert(s != NULL);

    smoother = s->smoother;
    A = s->A;
    b = s->b;
    x = s->x;
    nsweeps = s->nsweeps;
    istart = s->istart;
    iend = s->iend;
    istep = s->istep;
    relax = s->relax;
    ndeg = s->ndeg;
    order = s->cf_order;
    ordering = s->ordering;

    switch (smoother) {
        case SX_SM_GS:
            if (order && ordering != NULL) {
                sx_amg_smoother_gs_cf(x, A, b, nsweeps, ordering, -1);
            }
            else {
                sx_amg_smoother_gs(x, iend, istart, istep, A, b, nsweeps);
            }
            break;

        case SX_SM_SGS:
            sx_amg_smoother_sgs(x, A, b, nsweeps);
            break;

        case SX_SM_JACOBI:
            sx_amg_smoother_jacobi(x, iend, istart, istep, A, b, nsweeps);
            break;

        case SX_SM_L1DIAG:
            sx_amg_smoother_L1diag(x, iend, istart, istep, A, b, nsweeps);
            break;

        case SX_SM_POLY:
            sx_amg_smoother_poly(A, b, x, iend + 1, ndeg, nsweeps);
            break;

        case SX_SM_SOR:
            sx_amg_smoother_sor(x, iend, istart, istep, A, b, nsweeps, relax);
            break;

        case SX_SM_SSOR:
            sx_amg_smoother_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            sx_amg_smoother_sor(x, iend, istart, istep, A, b, nsweeps, relax);
            break;

        case SX_SM_GSOR:
            sx_amg_smoother_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            sx_amg_smoother_gs(x, iend, istart, istep, A, b, nsweeps);
            break;

        case SX_SM_SGSOR:
            sx_amg_smoother_sor(x, istart, iend, -istep, A, b, nsweeps, relax);
            sx_amg_smoother_sor(x, iend, istart, istep, A, b, nsweeps, relax);
            sx_amg_smoother_gs(x, istart, iend, -istep, A, b, nsweeps);
            sx_amg_smoother_gs(x, iend, istart, istep, A, b, nsweeps);
            break;

        default:
            sx_printf("### ERROR: Wrong smoother type %"dFMT"!\n", smoother);
            sx_exit_on_errcode(ERROR_INPUT_PAR, __FUNCTION__);
    }
}
