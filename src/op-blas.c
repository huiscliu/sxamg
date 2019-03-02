
#include "op-blas.h"
#include "internal.h"

/**
 * \fn void sx_blas_array_set (SX_INT n, SX_FLOAT *x, SX_FLOAT val)
 *
 * \brief Set initial value for an array to be x=val
 *
 * \pars n    Number of variables
 * \pars x    Pointer to the vector
 * \pars val  Initial value for the SX_FLOAT array
 *
 */
void sx_blas_array_set( SX_INT n, SX_FLOAT *x, SX_FLOAT val)
{
    SX_INT i;

    if (val == 0.0) {
        memset(x, 0x0, sizeof(SX_FLOAT) * n);
    }
    else {
        for (i = 0; i < n; ++i) x[i] = val;
    }
}

/**
 * \fn void sx_blas_array_cp (SX_INT n, const SX_FLOAT *x, SX_FLOAT *y) 
 *
 * \brief Copy an array to the other y=x
 *
 * \pars n    Number of variables
 * \pars x    Pointer to the original vector
 * \pars y    Pointer to the destination vector
 *
 */
void sx_blas_array_cp(SX_INT n, const SX_FLOAT *x, SX_FLOAT *y)
{
    memcpy(y, x, n * sizeof(SX_FLOAT));
}

void sx_blas_array_ax(SX_INT n, SX_FLOAT a, SX_FLOAT *x)
{
    SX_INT i;

    if (a == 1.0) {
        return;
    }
    else {
        for (i = 0; i < n; ++i) x[i] *= a;
    }
}

void sx_blas_array_axpy(SX_INT n, SX_FLOAT a, const SX_FLOAT *x, SX_FLOAT *y)
{
    SX_INT i;

    if (a == 1.0) {
        for (i = 0; i < n; ++i) y[i] += x[i];
    }
    else if (a == -1.0) {
        for (i = 0; i < n; ++i) y[i] -= x[i];
    }
    else {
        for (i = 0; i < n; ++i) y[i] += a * x[i];
    }
}

void sx_blas_array_axpyz(SX_INT n, SX_FLOAT a, const SX_FLOAT *x, const SX_FLOAT *y, SX_FLOAT *z)
{
    SX_INT i;

    for (i = 0; i < n; ++i) z[i] = a * x[i] + y[i];
}

void sx_blas_array_axpby(SX_INT n, SX_FLOAT a, const SX_FLOAT *x, SX_FLOAT b, SX_FLOAT *y)
{
    SX_INT i;

    for (i = 0; i < n; ++i) y[i] = a * x[i] + b * y[i];
}

SX_FLOAT sx_blas_array_dot(SX_INT n, const SX_FLOAT *x, const SX_FLOAT *y)
{
    SX_INT i;
    SX_FLOAT value = 0.0;

    for (i = 0; i < n; ++i) value += x[i] * y[i];

    return value;
}

SX_FLOAT sx_blas_array_norm1(SX_INT n, const SX_FLOAT *x)
{
    SX_INT i;
    SX_FLOAT onenorm = 0.;

    for (i = 0; i < n; ++i) onenorm += SX_ABS(x[i]);

    return onenorm;
}

SX_FLOAT sx_blas_array_norm2(SX_INT n, const SX_FLOAT *x)
{
    SX_INT i;
    SX_FLOAT twonorm = 0.;

    for (i = 0; i < n; ++i) twonorm += x[i] * x[i];

    return sqrt(twonorm);
}

SX_FLOAT sx_blas_array_norminf(SX_INT n, const SX_FLOAT *x)
{
    SX_INT i;
    SX_FLOAT infnorm = 0.0;

    for (i = 0; i < n; ++i) infnorm = SX_MAX(infnorm, SX_ABS(x[i]));

    return infnorm;
}

/**
 * \fn void sx_blas_vec_axpy (SX_FLOAT a, const SX_VEC *x, SX_VEC *y)
 *
 * \brief y = a*x + y
 *
 * \pars a   SX_FLOAT factor a
 * \pars x   Pointer to SX_VEC x
 * \pars y   Pointer to SX_VEC y
 *
 */
void sx_blas_vec_axpy(SX_FLOAT a, const SX_VEC *x, SX_VEC *y)
{
    SX_INT i, m = x->n;
    SX_FLOAT *xpt = x->d, *ypt = y->d;

    if ((y->n - m) != 0) {
        sx_printf("### ERROR: Two vectors have different dimensions!\n");
        sx_exit_on_errcode(ERROR_DATA_STRUCTURE, __FUNCTION__);
    }

    for (i = 0; i < m; ++i) ypt[i] += a * xpt[i];
}

/**
 * \fn void sx_blas_vec_axpby (SX_FLOAT a, const SX_VEC *x, SX_FLOAT b, SX_VEC *y)
 *
 * \brief y = a * x + b * y
 *
 * \pars a   SX_FLOAT factor a
 * \pars x   Pointer to SX_VEC x
 * \pars b   SX_FLOAT factor b
 * \pars y   Pointer to SX_VEC y
 *
 */
void sx_blas_vec_axpby(SX_FLOAT a, const SX_VEC *x, SX_FLOAT b, SX_VEC *y)
{
    SX_INT i, m = x->n;
    SX_FLOAT *xpt = x->d, *ypt = y->d;

    if ((y->n - m) != 0) {
        sx_printf("### ERROR: Two vectors have different dimensions!\n");
        sx_exit_on_errcode(ERROR_DATA_STRUCTURE, __FUNCTION__);
    }

    for (i = 0; i < m; ++i) ypt[i] = a * xpt[i] + b * ypt[i];
}

/**
 * \fn void sx_blas_vec_axpyz (SX_FLOAT a, const SX_VEC *x, const SX_VEC *y, SX_VEC *z) 
 *
 * \brief z = a*x + y, z is a third vector (z is cleared)
 *
 * \pars a   SX_FLOAT factor a
 * \pars x   Pointer to SX_VEC x
 * \pars y   Pointer to SX_VEC y
 * \pars z   Pointer to SX_VEC z
 *
 */
void sx_blas_vec_axpyz(SX_FLOAT a, const SX_VEC *x, const SX_VEC *y, SX_VEC * z)
{
    const SX_INT m = x->n;
    SX_FLOAT *xpt = x->d, *ypt = y->d, *zpt = z->d;

    if ((y->n - m) != 0) {
        sx_printf("### ERROR: Two vectors have different dimensions!\n");
        sx_exit_on_errcode(ERROR_DATA_STRUCTURE, __FUNCTION__);
    }

    assert(z->n == m);

    memcpy(ypt, zpt, m * sizeof(SX_VEC));
    sx_blas_array_axpy(m, a, xpt, zpt);
}

/**
 * \fn void sx_blas_vec_axpybz (SX_FLOAT a, const SX_VEC *x, SX_FLOAT b, const SX_VEC *y, SX_VEC *z) 
 *
 * \brief z = a*x + b*y, z is a third vector (z is cleared)
 *
 * \pars a   SX_FLOAT factor a
 * \pars x   Pointer to SX_VEC x
 * \pars b   SX_FLOAT factor b
 * \pars y   Pointer to SX_VEC y
 * \pars z   Pointer to SX_VEC z
 *
 */
void sx_blas_vec_axpbyz(SX_FLOAT a, const SX_VEC *x, SX_FLOAT b, const SX_VEC *y, SX_VEC * z)
{
    const SX_INT m = x->n;
    SX_INT i;
    SX_FLOAT *xpt = x->d, *ypt = y->d, *zpt = z->d;

    if ((y->n - m) != 0) {
        sx_printf("### ERROR: Two vectors have different dimensions!\n");
        sx_exit_on_errcode(ERROR_DATA_STRUCTURE, __FUNCTION__);
    }

    assert(z->n == m);

    for (i = 0; i < m; ++i) zpt[i] = a * xpt[i] + b * ypt[i];
}

/**
 * \fn SX_FLOAT sx_blas_vec_dot(const SX_VEC *x, const SX_VEC *y) 
 *
 * \brief Inner product of two vectors (x,y)
 *
 * \pars x   Pointer to SX_VEC x
 * \pars y   Pointer to SX_VEC y
 *
 * \return Inner product
 *
 */

SX_FLOAT sx_blas_vec_dot(const SX_VEC *x, const SX_VEC *y)
{
    SX_FLOAT value = 0;
    SX_INT i;
    const SX_INT length = x->n;
    SX_FLOAT *xpt = x->d, *ypt = y->d;

    for (i = 0; i < length; ++i) value += xpt[i] * ypt[i];

    return value;
}

/**
 * \fn SX_FLOAT sx_blas_vec_norm2 (const SX_VEC *x) 
 *
 * \brief L2 norm of SX_VEC x
 *
 * \pars x   Pointer to SX_VEC x
 *
 * \return L2 norm of x
 *
 */

SX_FLOAT sx_blas_vec_norm2(const SX_VEC *x)
{
    SX_FLOAT twonorm = 0;
    SX_INT i;
    const SX_INT length = x->n;
    SX_FLOAT *xpt = x->d;

    for (i = 0; i < length; ++i) twonorm += xpt[i] * xpt[i];

    return sqrt(twonorm);
}

void sx_blas_mat_mxy(const SX_MAT *A, const SX_VEC *x, SX_VEC *y)
{
    const SX_INT m = A->num_rows;
    const SX_INT *ia = A->Ap, *ja = A->Aj;
    const SX_FLOAT *aj = A->Ax;
    SX_INT i, k, begin_row, end_row, nnz_num_row;
    SX_FLOAT temp;

    for (i = 0; i < m; ++i) {
        temp = 0.0;
        begin_row = ia[i];
        end_row = ia[i + 1];
        nnz_num_row = end_row - begin_row;
        switch (nnz_num_row) {
            case 3:
                k = begin_row;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                break;
            case 4:
                k = begin_row;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                break;
            case 5:
                k = begin_row;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                break;
            case 6:
                k = begin_row;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                break;
            case 7:
                k = begin_row;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                k++;
                temp += aj[k] * x->d[ja[k]];
                break;
            default:
                for (k = begin_row; k < end_row; ++k) {
                    temp += aj[k] * x->d[ja[k]];
                }
                break;
        }
        y->d[i] = temp;
    }
}

void sx_blas_mat_amxpy(SX_FLOAT alpha, const SX_MAT *A, const SX_VEC *x, SX_VEC *y)
{
    const SX_INT m = A->num_rows;
    const SX_INT *ia = A->Ap, *ja = A->Aj;
    const SX_FLOAT *aj = A->Ax;
    SX_INT i, k, begin_row, end_row;
    SX_FLOAT temp;

    if (alpha == 1.0) {
        for (i = 0; i < m; ++i) {
            temp = 0.0;
            begin_row = ia[i];
            end_row = ia[i + 1];

            for (k = begin_row; k < end_row; ++k) temp += aj[k] * x->d[ja[k]];

            y->d[i] += temp;
        }
    }
    else if (alpha == -1.0) {
        for (i = 0; i < m; ++i) {
            temp = 0.0;
            begin_row = ia[i];
            end_row = ia[i + 1];

            for (k = begin_row; k < end_row; ++k) temp += aj[k] * x->d[ja[k]];

            y->d[i] -= temp;
        }
    }
    else {
        for (i = 0; i < m; ++i) {
            temp = 0.0;
            begin_row = ia[i];
            end_row = ia[i + 1];

            for (k = begin_row; k < end_row; ++k) temp += aj[k] * x->d[ja[k]];

            y->d[i] += temp * alpha;
        }
    }
}

SX_MAT sx_blas_mat_rap(const SX_MAT *R, const SX_MAT *A, const SX_MAT *P)
{
    SX_INT n_coarse = R->num_rows;
    SX_INT *R_i = R->Ap;
    SX_INT *R_j = R->Aj;
    SX_FLOAT *R_data = R->Ax;

    SX_INT n_fine = A->num_rows;
    SX_INT *A_i = A->Ap;
    SX_INT *A_j = A->Aj;
    SX_FLOAT *A_data = A->Ax;

    SX_INT *P_i = P->Ap;
    SX_INT *P_j = P->Aj;
    SX_FLOAT *P_data = P->Ax;

    SX_INT RAP_size;
    SX_INT *RAP_p = NULL;
    SX_INT *RAP_j = NULL;
    SX_FLOAT *RAP_x = NULL;

    SX_INT *Ps_marker = NULL;
    SX_INT *As_marker = NULL;

    SX_INT ic, i1, i2, i3, jj1, jj2, jj3;
    SX_INT jj_cnter, jj_row_begining;
    SX_FLOAT r_entry, r_a_product, r_a_p_product;

    SX_INT coarse_mul_nthreads = n_coarse;
    SX_INT fine_mul_nthreads = n_fine;
    SX_INT coarse_add_nthreads = n_coarse + 1;
    SX_INT minus_one_length = coarse_mul_nthreads + fine_mul_nthreads;
    SX_INT total_calloc = minus_one_length + coarse_add_nthreads + 1;
    SX_MAT RAP;

    Ps_marker = (SX_INT *) sx_mem_calloc(total_calloc, sizeof(SX_INT));
    As_marker = Ps_marker + coarse_mul_nthreads;

    /*------------------------------------------------------*
     *  First Pass: Determine size of RAP and set up RAP_p  *
     *------------------------------------------------------*/
    RAP_p = (SX_INT *) sx_mem_calloc(n_coarse + 1, sizeof(SX_INT));

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

    RAP_j = (SX_INT *) sx_mem_calloc(RAP_size, sizeof(SX_INT));
    RAP_x = (SX_FLOAT *) sx_mem_calloc(RAP_size, sizeof(SX_FLOAT));

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

    sx_mem_free(Ps_marker);

    return RAP;
}
