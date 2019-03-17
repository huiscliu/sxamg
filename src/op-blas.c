
#include "op-blas.h"
#include "internal.h"

/**
 * \fn void sx_blas_array_set (SX_INT n, SX_FLT *x, SX_FLT val)
 *
 * \brief Set initial value for an array to be x=val
 *
 * \pars n    Number of variables
 * \pars x    Pointer to the vector
 * \pars val  Initial value for the SX_FLT array
 *
 */
void sx_blas_array_set( SX_INT n, SX_FLT *x, SX_FLT val)
{
    SX_INT i;

    for (i = 0; i < n; ++i) x[i] = val;
}

/**
 * \fn void sx_blas_array_cp (SX_INT n, const SX_FLT *x, SX_FLT *y) 
 *
 * \brief Copy an array to the other y=x
 *
 * \pars n    Number of variables
 * \pars x    Pointer to the original vector
 * \pars y    Pointer to the destination vector
 *
 */
void sx_blas_array_cp(SX_INT n, const SX_FLT *x, SX_FLT *y)
{
    memcpy(y, x, n * sizeof(SX_FLT));
}

void sx_blas_array_ax(SX_INT n, SX_FLT a, SX_FLT *x)
{
    SX_INT i;

    for (i = 0; i < n; ++i) x[i] *= a;
}

void sx_blas_array_axpy(SX_INT n, SX_FLT a, const SX_FLT *x, SX_FLT *y)
{
    SX_INT i;

    for (i = 0; i < n; ++i) y[i] += a * x[i];
}

void sx_blas_array_axpyz(SX_INT n, SX_FLT a, const SX_FLT *x, const SX_FLT *y, SX_FLT *z)
{
    SX_INT i;

    for (i = 0; i < n; ++i) z[i] = a * x[i] + y[i];
}

void sx_blas_array_axpby(SX_INT n, SX_FLT a, const SX_FLT *x, SX_FLT b, SX_FLT *y)
{
    SX_INT i;

    for (i = 0; i < n; ++i) y[i] = a * x[i] + b * y[i];
}

SX_FLT sx_blas_array_dot(SX_INT n, const SX_FLT *x, const SX_FLT *y)
{
    SX_INT i;
    SX_FLT value = 0.0;

    for (i = 0; i < n; ++i) value += x[i] * y[i];

    return value;
}

SX_FLT sx_blas_array_norm1(SX_INT n, const SX_FLT *x)
{
    SX_INT i;
    SX_FLT onenorm = 0.;

    for (i = 0; i < n; ++i) onenorm += SX_ABS(x[i]);

    return onenorm;
}

SX_FLT sx_blas_array_norm2(SX_INT n, const SX_FLT *x)
{
    SX_INT i;
    SX_FLT twonorm = 0.;

    for (i = 0; i < n; ++i) twonorm += x[i] * x[i];

    return sqrt(twonorm);
}

SX_FLT sx_blas_array_norminf(SX_INT n, const SX_FLT *x)
{
    SX_INT i;
    SX_FLT infnorm = 0.0;

    for (i = 0; i < n; ++i) infnorm = SX_MAX(infnorm, SX_ABS(x[i]));

    return infnorm;
}

/**
 * \fn void sx_blas_vec_axpy (SX_FLT a, const SX_VEC *x, SX_VEC *y)
 *
 * \brief y = a*x + y
 *
 * \pars a   SX_FLT factor a
 * \pars x   Pointer to SX_VEC x
 * \pars y   Pointer to SX_VEC y
 *
 */
void sx_blas_vec_axpy(SX_FLT a, const SX_VEC *x, SX_VEC *y)
{
    SX_INT i, m = x->n;
    SX_FLT *xpt = x->d, *ypt = y->d;

    if ((y->n - m) != 0) {
        sx_printf("### ERROR: Two vectors have different dimensions!\n");
        sx_exit_on_errcode(ERROR_DATA_STRUCTURE, __FUNCTION__);
    }

    for (i = 0; i < m; ++i) ypt[i] += a * xpt[i];
}

/**
 * \fn void sx_blas_vec_axpby (SX_FLT a, const SX_VEC *x, SX_FLT b, SX_VEC *y)
 *
 * \brief y = a * x + b * y
 *
 * \pars a   SX_FLT factor a
 * \pars x   Pointer to SX_VEC x
 * \pars b   SX_FLT factor b
 * \pars y   Pointer to SX_VEC y
 *
 */
void sx_blas_vec_axpby(SX_FLT a, const SX_VEC *x, SX_FLT b, SX_VEC *y)
{
    SX_INT i, m = x->n;
    SX_FLT *xpt = x->d, *ypt = y->d;

    if ((y->n - m) != 0) {
        sx_printf("### ERROR: Two vectors have different dimensions!\n");
        sx_exit_on_errcode(ERROR_DATA_STRUCTURE, __FUNCTION__);
    }

    for (i = 0; i < m; ++i) ypt[i] = a * xpt[i] + b * ypt[i];
}

/**
 * \fn void sx_blas_vec_axpyz (SX_FLT a, const SX_VEC *x, const SX_VEC *y, SX_VEC *z) 
 *
 * \brief z = a*x + y, z is a third vector (z is cleared)
 *
 * \pars a   SX_FLT factor a
 * \pars x   Pointer to SX_VEC x
 * \pars y   Pointer to SX_VEC y
 * \pars z   Pointer to SX_VEC z
 *
 */
void sx_blas_vec_axpyz(SX_FLT a, const SX_VEC *x, const SX_VEC *y, SX_VEC * z)
{
    SX_INT m = x->n;
    SX_FLT *xpt = x->d, *ypt = y->d, *zpt = z->d;

    if ((y->n - m) != 0) {
        sx_printf("### ERROR: Two vectors have different dimensions!\n");
        sx_exit_on_errcode(ERROR_DATA_STRUCTURE, __FUNCTION__);
    }

    assert(z->n == m);

    memcpy(ypt, zpt, m * sizeof(SX_VEC));
    sx_blas_array_axpy(m, a, xpt, zpt);
}

/**
 * \fn void sx_blas_vec_axpybz (SX_FLT a, const SX_VEC *x, SX_FLT b, const SX_VEC *y, SX_VEC *z) 
 *
 * \brief z = a*x + b*y, z is a third vector (z is cleared)
 *
 * \pars a   SX_FLT factor a
 * \pars x   Pointer to SX_VEC x
 * \pars b   SX_FLT factor b
 * \pars y   Pointer to SX_VEC y
 * \pars z   Pointer to SX_VEC z
 *
 */
void sx_blas_vec_axpbyz(SX_FLT a, const SX_VEC *x, SX_FLT b, const SX_VEC *y, SX_VEC * z)
{
    SX_INT m = x->n;
    SX_INT i;
    SX_FLT *xpt = x->d, *ypt = y->d, *zpt = z->d;

    if ((y->n - m) != 0) {
        sx_printf("### ERROR: Two vectors have different dimensions!\n");
        sx_exit_on_errcode(ERROR_DATA_STRUCTURE, __FUNCTION__);
    }

    assert(z->n == m);

    for (i = 0; i < m; ++i) zpt[i] = a * xpt[i] + b * ypt[i];
}

/**
 * \fn SX_FLT sx_blas_vec_dot(const SX_VEC *x, const SX_VEC *y) 
 *
 * \brief Inner product of two vectors (x,y)
 *
 * \pars x   Pointer to SX_VEC x
 * \pars y   Pointer to SX_VEC y
 *
 * \return Inner product
 *
 */

SX_FLT sx_blas_vec_dot(const SX_VEC *x, const SX_VEC *y)
{
    SX_FLT value = 0;
    SX_INT i;
    SX_INT length = x->n;
    SX_FLT *xpt = x->d, *ypt = y->d;

    for (i = 0; i < length; ++i) value += xpt[i] * ypt[i];

    return value;
}

/**
 * \fn SX_FLT sx_blas_vec_norm2 (const SX_VEC *x) 
 *
 * \brief L2 norm of SX_VEC x
 *
 * \pars x   Pointer to SX_VEC x
 *
 * \return L2 norm of x
 *
 */

SX_FLT sx_blas_vec_norm2(const SX_VEC *x)
{
    SX_FLT twonorm = 0;
    SX_INT i;
    SX_INT length = x->n;
    SX_FLT *xpt = x->d;

    for (i = 0; i < length; ++i) twonorm += xpt[i] * xpt[i];

    return sqrt(twonorm);
}

void sx_blas_vec_copy(const SX_VEC *x, SX_VEC *y)
{
    assert(x != NULL);
    assert(y != NULL);
    assert(x->n == y->n);
    assert(x->n >= 0);

    memcpy(y->d, x->d, x->n * sizeof(*x->d));
}

void sx_blas_vec_set(SX_VEC *x, SX_FLT val)
{
    int i;

    assert(x != NULL);
    assert(x->n >= 0);

    for (i = 0; i < x->n; i++) x->d[i] = val;
}

void sx_blas_mv_mxy(const SX_MAT *A, const SX_VEC *x, SX_VEC *y)
{
    SX_INT m = A->num_rows;
    SX_INT *ia = A->Ap, *ja = A->Aj;
    SX_FLT *aj = A->Ax;
    SX_INT i, k, begin_row, end_row;
    SX_FLT temp;

    for (i = 0; i < m; ++i) {
        temp = 0.0;
        begin_row = ia[i];
        end_row = ia[i + 1];

        for (k = begin_row; k < end_row; ++k) {
            temp += aj[k] * x->d[ja[k]];
        }

        y->d[i] = temp;
    }
}

void sx_blas_mv_amxpy(SX_FLT alpha, const SX_MAT *A, const SX_VEC *x, SX_VEC *y)
{
    SX_INT m = A->num_rows;
    SX_INT *ia = A->Ap, *ja = A->Aj;
    SX_FLT *aj = A->Ax;
    SX_INT i, k, begin_row, end_row;
    SX_FLT temp;

    for (i = 0; i < m; ++i) {
        temp = 0.0;
        begin_row = ia[i];
        end_row = ia[i + 1];

        for (k = begin_row; k < end_row; ++k) temp += aj[k] * x->d[ja[k]];

        y->d[i] += temp * alpha;
    }
}

void sx_blas_mv_amxpby(SX_FLT alpha, const SX_MAT *A, const SX_VEC *x, SX_FLT beta, SX_VEC *y)
{
    SX_INT m = A->num_rows;
    SX_INT *ia = A->Ap, *ja = A->Aj;
    SX_FLT *aj = A->Ax;
    SX_INT i, k, begin_row, end_row;
    SX_FLT temp;

    for (i = 0; i < m; ++i) {
        temp = 0.0;
        begin_row = ia[i];
        end_row = ia[i + 1];

        for (k = begin_row; k < end_row; ++k) temp += aj[k] * x->d[ja[k]];

        y->d[i] = temp * alpha + y->d[i] * beta;
    }
}

void sx_blas_mv_amxpbyz(SX_FLT alpha, const SX_MAT *A, const SX_VEC *x, SX_FLT beta, SX_VEC *y, SX_VEC *z)
{
    SX_INT m = A->num_rows;
    SX_INT *ia = A->Ap, *ja = A->Aj;
    SX_FLT *aj = A->Ax;
    SX_INT i, k, begin_row, end_row;
    SX_FLT temp;

    for (i = 0; i < m; ++i) {
        temp = 0.0;
        begin_row = ia[i];
        end_row = ia[i + 1];

        for (k = begin_row; k < end_row; ++k) temp += aj[k] * x->d[ja[k]];

        z->d[i] = temp * alpha + y->d[i] * beta;
    }
}

SX_MAT sx_blas_mat_rap(const SX_MAT *R, const SX_MAT *A, const SX_MAT *P)
{
    SX_INT n_coarse = R->num_rows;
    SX_INT *R_i = R->Ap;
    SX_INT *R_j = R->Aj;
    SX_FLT *R_data = R->Ax;

    SX_INT n_fine = A->num_rows;
    SX_INT *A_i = A->Ap;
    SX_INT *A_j = A->Aj;
    SX_FLT *A_data = A->Ax;

    SX_INT *P_i = P->Ap;
    SX_INT *P_j = P->Aj;
    SX_FLT *P_data = P->Ax;

    SX_INT RAP_size;
    SX_INT *RAP_p = NULL;
    SX_INT *RAP_j = NULL;
    SX_FLT *RAP_x = NULL;

    SX_INT *Ps_marker = NULL;
    SX_INT *As_marker = NULL;

    SX_INT ic, i1, i2, i3, jj1, jj2, jj3;
    SX_INT jj_cnter, jj_row_begining;
    SX_FLT r_entry, r_a_product, r_a_p_product;

    SX_INT coarse_mul = n_coarse;
    SX_INT fine_mul = n_fine;
    SX_INT coarse_add = n_coarse + 1;
    SX_INT minus_one_length = coarse_mul + fine_mul;
    SX_INT total_calloc = minus_one_length + coarse_add + 1;
    SX_MAT RAP;

    Ps_marker = (SX_INT *) sx_calloc(total_calloc, sizeof(SX_INT));
    As_marker = Ps_marker + coarse_mul;

    /*------------------------------------------------------*
     *  First Pass: Determine size of RAP and set up RAP_p  *
     *------------------------------------------------------*/
    RAP_p = (SX_INT *) sx_calloc(n_coarse + 1, sizeof(SX_INT));

    sx_iarray_set(minus_one_length, Ps_marker, -1);

    jj_cnter = 0;
    for (ic = 0; ic < n_coarse; ic++) {
        Ps_marker[ic] = jj_cnter;
        jj_row_begining = jj_cnter;
        jj_cnter++;

        for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {
            i1 = R_j[jj1];

            for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {
                i2 = A_j[jj2];
                if (As_marker[i2] != ic) {
                    As_marker[i2] = ic;
                    for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                        i3 = P_j[jj3];
                        if (Ps_marker[i3] < jj_row_begining) {
                            Ps_marker[i3] = jj_cnter;
                            jj_cnter++;
                        }
                    }
                }
            }
        }

        RAP_p[ic] = jj_row_begining;
    }

    RAP_p[n_coarse] = jj_cnter;
    RAP_size = jj_cnter;

    RAP_j = (SX_INT *) sx_calloc(RAP_size, sizeof(SX_INT));
    RAP_x = (SX_FLT *) sx_calloc(RAP_size, sizeof(SX_FLT));

    sx_iarray_set(minus_one_length, Ps_marker, -1);

    jj_cnter = 0;
    for (ic = 0; ic < n_coarse; ic++) {
        Ps_marker[ic] = jj_cnter;
        jj_row_begining = jj_cnter;
        RAP_j[jj_cnter] = ic;
        RAP_x[jj_cnter] = 0.0;
        jj_cnter++;

        for (jj1 = R_i[ic]; jj1 < R_i[ic + 1]; jj1++) {
            r_entry = R_data[jj1];

            i1 = R_j[jj1];
            for (jj2 = A_i[i1]; jj2 < A_i[i1 + 1]; jj2++) {
                r_a_product = r_entry * A_data[jj2];

                i2 = A_j[jj2];
                if (As_marker[i2] != ic) {
                    As_marker[i2] = ic;

                    for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                        r_a_p_product = r_a_product * P_data[jj3];

                        i3 = P_j[jj3];
                        if (Ps_marker[i3] < jj_row_begining) {
                            Ps_marker[i3] = jj_cnter;
                            RAP_x[jj_cnter] = r_a_p_product;
                            RAP_j[jj_cnter] = i3;
                            jj_cnter++;
                        }
                        else {
                            RAP_x[Ps_marker[i3]] += r_a_p_product;
                        }
                    }
                }
                else {
                    for (jj3 = P_i[i2]; jj3 < P_i[i2 + 1]; jj3++) {
                        i3 = P_j[jj3];

                        r_a_p_product = r_a_product * P_data[jj3];
                        RAP_x[Ps_marker[i3]] += r_a_p_product;
                    }
                }
            }
        }
    }

    RAP.num_rows = n_coarse;
    RAP.num_cols = n_coarse;
    RAP.num_nnzs = RAP_size;
    RAP.Ap = RAP_p;
    RAP.Aj = RAP_j;
    RAP.Ax = RAP_x;

    sx_free(Ps_marker);

    return RAP;
}

/**
 * \fn void sx_blas_mat_mxm (const SX_MAT  *A, const SX_MAT  *B, SX_MAT *C)
 *
 * \brief Sparse matrix multiplication C=A*B
 *
 * \param A   Pointer to the dCSRmat matrix A
 * \param B   Pointer to the dCSRmat matrix B
 * \param C   Pointer to dCSRmat matrix equal to A*B
 *
 */
void sx_blas_mat_mxm(const SX_MAT * A, const SX_MAT * B, SX_MAT * C)
{
    SX_INT i, j, k, l, count;
    SX_INT countJD;
    SX_INT *JD = (SX_INT *) sx_calloc(B->num_cols, sizeof(SX_INT));

    C->num_rows = A->num_rows;
    C->num_cols = B->num_cols;
    C->Ax = NULL;
    C->Aj = NULL;
    C->Ap = (SX_INT *) sx_calloc(C->num_rows + 1, sizeof(SX_INT));

    for (i = 0; i < B->num_cols; ++i)
        JD[i] = -1;

    // step 1: Find first the structure IA of C
    for (i = 0; i < C->num_rows; ++i) {
        count = 0;

        for (k = A->Ap[i]; k < A->Ap[i + 1]; ++k) {
            for (j = B->Ap[A->Aj[k]]; j < B->Ap[A->Aj[k] + 1]; ++j) {
                for (l = 0; l < count; l++) {
                    if (JD[l] == B->Aj[j])
                        break;
                }

                if (l == count) {
                    JD[count] = B->Aj[j];
                    count++;
                }
            }
        }
        C->Ap[i + 1] = count;
        for (j = 0; j < count; ++j) {
            JD[j] = -1;
        }
    }

    for (i = 0; i < C->num_rows; ++i)
        C->Ap[i + 1] += C->Ap[i];

    // step 2: Find the structure JA of C
    C->Aj = (SX_INT *) sx_calloc(C->Ap[C->num_rows], sizeof(SX_INT));

    for (i = 0; i < C->num_rows; ++i) {
        countJD = 0;
        count = C->Ap[i];
        for (k = A->Ap[i]; k < A->Ap[i + 1]; ++k) {
            for (j = B->Ap[A->Aj[k]]; j < B->Ap[A->Aj[k] + 1]; ++j) {
                for (l = 0; l < countJD; l++) {
                    if (JD[l] == B->Aj[j])
                        break;
                }

                if (l == countJD) {
                    C->Aj[count] = B->Aj[j];
                    JD[countJD] = B->Aj[j];
                    count++;
                    countJD++;
                }
            }
        }

        //for (j=0;j<countJD;++j) JD[j]=-1;
        sx_iarray_set(countJD, JD, -1);
    }

    sx_free(JD);
    JD = NULL;

    // step 3: Find the structure A of C
    C->Ax = (SX_FLT *) sx_calloc(C->Ap[C->num_rows], sizeof(SX_FLT));

    for (i = 0; i < C->num_rows; ++i) {
        for (j = C->Ap[i]; j < C->Ap[i + 1]; ++j) {
            C->Ax[j] = 0;
            for (k = A->Ap[i]; k < A->Ap[i + 1]; ++k) {
                for (l = B->Ap[A->Aj[k]]; l < B->Ap[A->Aj[k] + 1]; l++) {
                    if (B->Aj[l] == C->Aj[j]) {
                        C->Ax[j] += A->Ax[k] * B->Ax[l];
                    }           // end if
                }               // end for l
            }                   // end for k
        }                       // end for j
    }                           // end for i

    C->num_nnzs = C->Ap[C->num_rows] - C->Ap[0];
}
