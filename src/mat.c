
#include "mat.h"
#include "internal.h"

extern void sx_iarray_cp(const SX_INT n,SX_INT *x,SX_INT *y);
extern void sx_iarray_set(const SX_INT n,SX_INT *x,const SX_INT val);

extern void sx_blas_array_cp(const SX_INT n,SX_FLT *x,SX_FLT *y);
extern void sx_blas_array_set(const SX_INT n,SX_FLT *x,const SX_FLT val);

/**
 * \fn SX_MAT sx_mat_struct_create (SX_INT m, SX_INT n, SX_INT nnz)
 *
 * \brief Create CSR sparse matrix data memory space
 *
 * \pars m    Number of rows
 * \pars n    Number of columns
 * \pars nnz  Number of nonzeros
 *
 * \return A   the new SX_MAT matrix
 *
 */
SX_MAT sx_mat_struct_create(SX_INT m, SX_INT n, SX_INT nnz)
{
    SX_MAT A;

    assert(m > 0);
    assert(n > 0);
    assert(nnz >= 0);

    A.Ap = (SX_INT *) sx_calloc(m + 1, sizeof(SX_INT));

    if (nnz > 0) {
        A.Aj = (SX_INT *) sx_calloc(nnz, sizeof(SX_INT));
    }
    else {
        A.Aj = NULL;
    }

    if (nnz > 0) {
        A.Ax = (SX_FLT *) sx_calloc(nnz, sizeof(SX_FLT));
    }
    else {
        A.Ax = NULL;
    }

    A.num_rows = m;
    A.num_cols = n;
    A.num_nnzs = nnz;

    return A;
}

/**
 * \fn SX_IMAT sx_imat_struct_create (SX_INT m, SX_INT n, SX_INT nnz)
 *
 * \brief Create CSR sparse matrix data memory space
 *
 * \pars m    Number of rows
 * \pars n    Number of columns
 * \pars nnz  Number of nonzeros
 *
 * \return A   the new SX_IMAT matrix
 *
 */
SX_IMAT sx_imat_struct_create(SX_INT m, SX_INT n, SX_INT nnz)
{
    SX_IMAT A;

    assert(m > 0);
    assert(n > 0);
    assert(nnz >= 0);

    A.Ap = (SX_INT *) sx_calloc(m + 1, sizeof(SX_INT));

    if (nnz > 0) {
        A.Aj = (SX_INT *) sx_calloc(nnz, sizeof(SX_INT));
    }
    else {
        A.Aj = NULL;
    }

    if (nnz > 0) {
        A.Ax = (SX_INT *) sx_calloc(nnz, sizeof(SX_INT));
    }
    else {
        A.Ax = NULL;
    }

    A.num_rows = m;
    A.num_cols = n;
    A.num_nnzs = nnz;

    return A;
}

SX_MAT sx_mat_create(SX_INT nrow, SX_INT ncol, SX_INT *Ap, SX_INT *Aj, SX_FLT *Ax)
{
    SX_MAT A;
    SX_INT nnz;

    assert(nrow > 0);
    assert(Ap != NULL);

    nnz = Ap[nrow];
    A = sx_mat_struct_create(nrow, ncol, nnz);

    if (nnz > 0) {
        assert(Aj != NULL);
        assert(Ax != NULL);
    }

    memcpy(A.Ap, Ap, (nrow + 1) * sizeof(SX_INT));
    memcpy(A.Aj, Aj, nnz * sizeof(SX_INT));
    memcpy(A.Ax, Ax, nnz * sizeof(SX_FLT));

    return A;
}

/**
 * \fn void sx_mat_destroy (SX_MAT *A)
 *
 * \brief Free CSR sparse matrix data memory space
 *
 * \pars A   Pointer to the SX_MAT matrix
 *
 */
void sx_mat_destroy(SX_MAT *A)
{
    if (A == NULL) return;

    sx_free(A->Ap);
    A->Ap = NULL;

    sx_free(A->Aj);
    A->Aj = NULL;

    sx_free(A->Ax);
    A->Ax = NULL;
}

/**
 * \fn void sx_imat_destroy (SX_IMAT *A)
 *
 * \brief Free CSR sparse matrix data memory space
 *
 * \pars A   Pointer to the SX_IMAT matrix
 *
 */
void sx_imat_destroy(SX_IMAT *A)
{
    if (A == NULL) return;

    sx_free(A->Ap);
    A->Ap = NULL;

    sx_free(A->Aj);
    A->Aj = NULL;

    sx_free(A->Ax);
    A->Ax = NULL;
}

/**
 * \fn static void iSwapping (SX_INT *w, const SX_INT i, const SX_INT j)
 *
 * \brief swap the i-th and j-th element in the array 'w' (SX_INT type)
 *
 * \pars w    Pointer to the array
 * \pars i    One entry in w
 * \pars j    Another entry in w
 *
 */
static void iSwapping(SX_INT *w, const SX_INT i, const SX_INT j)
{
    SX_INT temp = w[i];
    w[i] = w[j];
    w[j] = temp;
}

/*!
 * \fn void sx_aux_iQuickSortIndex (SX_INT *a, SX_INT left, SX_INT right, SX_INT *index)
 *
 * \brief Reorder the index of (SX_INT type) so that 'a' is in ascending order
 *
 * \pars a       Pointer to the array
 * \pars left    Starting index
 * \pars right   Ending index
 * \pars index   Index of 'a' (out)
 *
 * \note 'left' and 'right' are usually set to be 0 and n-1,respectively,where n is the
 *       length of 'a'. 'index' should be initialized in the nature order and it has the
 *       same length as 'a'.
 */
static void sx_aux_iQuickSortIndex(SX_INT *a, SX_INT left, SX_INT right, SX_INT *index)
{
    SX_INT i, last;

    if (left >= right)
        return;

    iSwapping(index, left, (left + right) / 2);

    last = left;
    for (i = left + 1; i <= right; ++i) {
        if (a[index[i]] < a[index[left]]) {
            iSwapping(index, ++last, i);
        }
    }

    iSwapping(index, left, last);

    sx_aux_iQuickSortIndex(a, left, last - 1, index);
    sx_aux_iQuickSortIndex(a, last + 1, right, index);
}

/**
 * \fn void sx_mat_sort (SX_MAT *A)
 *
 * \brief Sort each row of A in ascending order w.r.t. column indices.
 *
 * \pars A   Pointer to the SX_MAT matrix
 *
 */
void sx_mat_sort(SX_MAT *A)
{
    const SX_INT n = A->num_cols;
    SX_INT i, j, start, row_length;

    // temp memory for sorting rows of A
    SX_INT *index, *ja;
    SX_FLT *a;

    index = (SX_INT *) sx_calloc(n, sizeof(SX_INT));
    ja = (SX_INT *) sx_calloc(n, sizeof(SX_INT));
    a = (SX_FLT *) sx_calloc(n, sizeof(SX_FLT));

    for (i = 0; i < n; ++i) {
        start = A->Ap[i];
        row_length = A->Ap[i + 1] - start;

        for (j = 0; j < row_length; ++j)
            index[j] = j;

        sx_aux_iQuickSortIndex(&(A->Aj[start]), 0, row_length - 1, index);

        for (j = 0; j < row_length; ++j) {
            ja[j] = A->Aj[start + index[j]];
            a[j] = A->Ax[start + index[j]];
        }

        for (j = 0; j < row_length; ++j) {
            A->Aj[start + j] = ja[j];
            A->Ax[start + j] = a[j];
        }
    }

    // clean up memory
    sx_free(index);
    sx_free(ja);
    sx_free(a);
}

/**
 * \fn SX_VEC sx_mat_get_diag (SX_MAT *A, SX_INT n)
 *
 * \brief Get first n diagonal entries of a CSR matrix A
 *
 * \pars A     Pointer to SX_MAT CSR matrix
 * \pars n     Number of diagonal entries to get (if n=0, then get all diagonal entries)
 *
 */
SX_VEC sx_mat_get_diag(SX_MAT *A, SX_INT n)
{
    SX_INT i, k, j, ibegin, iend;
    SX_VEC diag;

    if (n == 0 || n > A->num_rows || n > A->num_cols)
        n = SX_MIN(A->num_rows, A->num_cols);

    diag.n = n;
    diag.d = (SX_FLT *) sx_calloc(n, sizeof(SX_FLT));

    for (i = 0; i < n; ++i) {
        ibegin = A->Ap[i];
        iend = A->Ap[i + 1];
        for (k = ibegin; k < iend; ++k) {
            j = A->Aj[k];
            if ((j - i) == 0) {
                diag.d[i] = A->Ax[k];
                break;
            }
        }
    }

    return diag;
}

/**
 * \fn void sx_imat_cp (SX_IMAT *src, SX_IMAT *des)
 *
 * \brief Copy a SX_IMAT to a new one des=src
 *
 * \pars src   Pointer to the SX_IMAT matrix
 * \pars des   Pointer to the SX_IMAT matrix
 *
 */
void sx_imat_cp(SX_IMAT *src, SX_IMAT *des)
{
    des->num_rows = src->num_rows;
    des->num_cols = src->num_cols;
    des->num_nnzs = src->num_nnzs;

    sx_iarray_cp(src->num_rows + 1, src->Ap, des->Ap);
    sx_iarray_cp(src->num_nnzs, src->Aj, des->Aj);
    sx_iarray_cp(src->num_nnzs, src->Ax, des->Ax);
}

/**
 * \fn void sx_mat_cp (SX_MAT *src, SX_MAT *des)
 *
 * \brief copy a SX_MAT to a new one des=src
 *
 * \pars src   Pointer to the SX_MAT matrix
 * \pars des   Pointer to the SX_MAT matrix
 *
 */
void sx_mat_cp(SX_MAT *src, SX_MAT *des)
{
    des->num_rows = src->num_rows;
    des->num_cols = src->num_cols;
    des->num_nnzs = src->num_nnzs;

    sx_iarray_cp(src->num_rows + 1, src->Ap, des->Ap);
    sx_iarray_cp(src->num_nnzs, src->Aj, des->Aj);
    sx_blas_array_cp(src->num_nnzs, src->Ax, des->Ax);
}

/**
 * \fn SX_IMAT sx_imat_trans (SX_IMAT *A)
 *
 * \brief Find transpose of SX_IMAT matrix A
 *
 * \pars A   Pointer to the SX_IMAT matrix A
 *
 * \return    The transpose of SX_IMAT matrix A
 *
 */
SX_IMAT sx_imat_trans(SX_IMAT *A)
{
    const SX_INT n = A->num_rows, m = A->num_cols, nnz = A->num_nnzs, m1 = m - 1;

    // Local variables
    SX_INT i, j, k, p;
    SX_INT ibegin, iend;
    SX_IMAT AT;

    AT.num_rows = m;
    AT.num_cols = n;
    AT.num_nnzs = nnz;

    AT.Ap = (SX_INT *) sx_calloc(m + 1, sizeof(SX_INT));

    AT.Aj = (SX_INT *) sx_calloc(nnz, sizeof(SX_INT));

    if (A->Ax) {
        AT.Ax = (SX_INT *) sx_calloc(nnz, sizeof(SX_INT));
    }
    else {
        AT.Ax = NULL;
    }

    // first pass: find the Number of nonzeros in the first m-1 columns of A
    // Note: these Numbers are stored in the array AT.Ap from 1 to m-1
    sx_iarray_set(m + 1, AT.Ap, 0);

    for (j = 0; j < nnz; ++j) {
        i = A->Aj[j];           // column Number of A = row Number of A'
        if (i < m1) AT.Ap[i + 2]++;
    }

    for (i = 2; i <= m; ++i) AT.Ap[i] += AT.Ap[i - 1];

    // second pass: form A'
    if (A->Ax != NULL) {
        for (i = 0; i < n; ++i) {
            ibegin = A->Ap[i], iend = A->Ap[i + 1];

            for (p = ibegin; p < iend; p++) {
                j = A->Aj[p] + 1;
                k = AT.Ap[j];
                AT.Aj[k] = i;
                AT.Ax[k] = A->Ax[p];
                AT.Ap[j] = k + 1;
            }
        }
    }
    else {
        for (i = 0; i < n; ++i) {
            ibegin = A->Ap[i], iend = A->Ap[i + 1];
            for (p = ibegin; p < iend; p++) {
                j = A->Aj[p] + 1;
                k = AT.Ap[j];
                AT.Aj[k] = i;
                AT.Ap[j] = k + 1;
            }
        }
    }

    return AT;
}

/**
 * \fn SX_MAT sx_mat_trans (SX_MAT *A)
 *
 * \brief Find transpose of SX_MAT matrix A
 *
 * \pars A   Pointer to the SX_MAT matrix
 *
 */
SX_MAT sx_mat_trans(SX_MAT *A)
{
    const SX_INT n = A->num_rows, m = A->num_cols, nnz = A->num_nnzs;
    SX_INT i, j, k, p;
    SX_MAT AT;

    AT.num_rows = m;
    AT.num_cols = n;
    AT.num_nnzs = nnz;

    AT.Ap = (SX_INT *) sx_calloc(m + 1, sizeof(SX_INT));
    AT.Aj = (SX_INT *) sx_calloc(nnz, sizeof(SX_INT));

    if (A->Ax) {
        AT.Ax = (SX_FLT *) sx_calloc(nnz, sizeof(SX_FLT));
    }
    else {
        AT.Ax = NULL;
    }

    // sx_iarray_set(m+1, AT.Ap, 0);
    memset(AT.Ap, 0, sizeof(SX_INT) * (m + 1));

    for (j = 0; j < nnz; ++j) {
        i = A->Aj[j];
        if (i < m - 1) AT.Ap[i + 2]++;
    }

    for (i = 2; i <= m; ++i) AT.Ap[i] += AT.Ap[i - 1];

    // second pass: form A'
    if (A->Ax) {
        for (i = 0; i < n; ++i) {
            SX_INT ibegin = A->Ap[i], iend = A->Ap[i + 1];
            for (p = ibegin; p < iend; p++) {
                j = A->Aj[p] + 1;
                k = AT.Ap[j];
                AT.Aj[k] = i;
                AT.Ax[k] = A->Ax[p];
                AT.Ap[j] = k + 1;
            }
        }
    }
    else {
        for (i = 0; i < n; ++i) {
            SX_INT ibegin = A->Ap[i], iend1 = A->Ap[i + 1];
            for (p = ibegin; p < iend1; p++) {
                j = A->Aj[p] + 1;
                k = AT.Ap[j];
                AT.Aj[k] = i;
                AT.Ap[j] = k + 1;
            }
        }
    }

    return AT;
}
