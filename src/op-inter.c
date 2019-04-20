
#include "op-inter.h"
#include "internal.h"

static void interp_DIR(SX_MAT *, SX_IVEC *, SX_MAT *, SX_AMG_PARS *);
static void interp_STD(SX_MAT *, SX_IVEC *, SX_MAT *, SX_IMAT *, SX_AMG_PARS *);

/**
 * \fn void sx_amg_interp(SX_MAT *A, SX_IVEC *vertices, SX_MAT *P, SX_IMAT *S,
 *                           SX_AMG_PARS *pars)
 *
 * \brief Generate interpolation operator P
 *
 * \pars A          Pointer to SX_MAT: the coefficient matrix (index starts from 0)
 * \pars vertices   Indicator vector for the C/F splitting of the variables
 * \pars P          Prolongation (input: nonzero pattern, output: prolongation)
 * \pars S          Strong connection matrix
 * \pars pars      AMG parseters
 *
 */
void sx_amg_interp(SX_MAT *A, SX_IVEC *vertices, SX_MAT *P, SX_IMAT *S, SX_AMG_PARS *pars)
{
    SX_INTERP_TYPE interp_type = pars->interp_type;

    switch (interp_type) {
        case SX_INTERP_DIR:       // Direct interpolation
            interp_DIR(A, vertices, P, pars);
            break;

        case SX_INTERP_STD:       // Standard interpolation
            interp_STD(A, vertices, P, S, pars);
            break;

        default:
            sx_exit_on_errcode(ERROR_AMG_INTERP_TYPE, __FUNCTION__);
    }
}

/**
 * \fn void sx_amg_interp_trunc(SX_MAT *P, SX_AMG_PARS *pars)
 *
 * \brief Truncation step for prolongation operators
 *
 * \pars P       Prolongation (input: full, output: truncated)
 * \pars pars    Pointer to SX_AMG_PARS: AMG parseters
 *
 */
void sx_amg_interp_trunc(SX_MAT *P, SX_AMG_PARS *pars)
{
    const SX_INT row = P->num_rows;
    const SX_INT nnzold = P->num_nnzs;
    const SX_INT verb = pars->verb;
    const SX_FLT eps_tr = pars->trunc_threshold;

    // local variables
    SX_INT num_nonzero = 0;        // number of non zeros after truncation
    SX_FLT Min_neg, Max_pos;     // min negative and max positive entries
    SX_FLT Fac_neg, Fac_pos;     // factors for negative and positive entries
    SX_FLT Sum_neg, TSum_neg;    // sum and truncated sum of negative entries
    SX_FLT Sum_pos, TSum_pos;    // sum and truncated sum of positive entries

    SX_INT index1 = 0, index2 = 0, begin, end;
    SX_INT i, j;

    for (i = 0; i < row; ++i) {
        begin = P->Ap[i];
        end = P->Ap[i + 1];

        P->Ap[i] = num_nonzero;
        Min_neg = Max_pos = 0;
        Sum_neg = Sum_pos = 0;
        TSum_neg = TSum_pos = 0;

        // 1. Summations of positive and negative entries
        for (j = begin; j < end; ++j) {
            if (P->Ax[j] > 0) {
                Sum_pos += P->Ax[j];
                Max_pos = SX_MAX(Max_pos, P->Ax[j]);
            }
            else if (P->Ax[j] < 0) {
                Sum_neg += P->Ax[j];
                Min_neg = SX_MIN(Min_neg, P->Ax[j]);
            }
        }

        Max_pos *= eps_tr;
        Min_neg *= eps_tr;

        // 2. Set JA of truncated P
        for (j = begin; j < end; ++j) {
            if (P->Ax[j] >= Max_pos) {
                num_nonzero++;
                P->Aj[index1++] = P->Aj[j];
                TSum_pos += P->Ax[j];
            }
            else if (P->Ax[j] <= Min_neg) {
                num_nonzero++;
                P->Aj[index1++] = P->Aj[j];
                TSum_neg += P->Ax[j];
            }
        }

        // 3. Compute factors and set values of truncated P
        if (TSum_pos > SMALLFLOAT) {
            Fac_pos = Sum_pos / TSum_pos;       // factor for positive entries
        }
        else {
            Fac_pos = 1.0;
        }

        if (TSum_neg < -SMALLFLOAT) {
            Fac_neg = Sum_neg / TSum_neg;       // factor for negative entries
        }
        else {
            Fac_neg = 1.0;
        }

        for (j = begin; j < end; ++j) {
            if (P->Ax[j] >= Max_pos)
                P->Ax[index2++] = P->Ax[j] * Fac_pos;
            else if (P->Ax[j] <= Min_neg)
                P->Ax[index2++] = P->Ax[j] * Fac_neg;
        }
    }

    // resize the truncated prolongation P
    P->num_nnzs = P->Ap[row] = num_nonzero;
    P->Aj = (SX_INT *) sx_realloc(P->Aj, num_nonzero * sizeof(SX_INT));
    P->Ax = (SX_FLT *) sx_realloc(P->Ax, num_nonzero * sizeof(SX_FLT));

    if (verb >= 4) {
        sx_printf("Truncate prolongation, nnz before: %10"dFMT", after: %10"dFMT"\n",
               nnzold, num_nonzero);
    }
}

/**
 * \fn static void interp_DIR(SX_MAT *A, SX_IVEC *vertices, SX_MAT *P, SX_AMG_PARS *pars)
 *
 * \brief Direct interpolation
 *
 * \pars A          Pointer to SX_MAT: the coefficient matrix (index starts from 0)
 * \pars vertices   Indicator vector for the C/F splitting of the variables
 * \pars P          Prolongation (input: nonzero pattern, output: prolongation)
 * \pars pars      Pointer to SX_AMG_PARS: AMG parseters
 *
 */
static void interp_DIR(SX_MAT * A, SX_IVEC * vertices, SX_MAT * P, SX_AMG_PARS * pars)
{
    SX_INT row = A->num_rows;
    SX_INT *vec = vertices->d;

    // local variables
    SX_INT IS_STRONG;              // is the variable strong coupled to i?
    SX_INT num_pcouple;            // number of positive strong couplings
    SX_INT begin_row, end_row;
    SX_INT i, j, k, l, index = 0, idiag;

    // a_minus and a_plus for Neighbors and Prolongation support
    SX_FLT amN, amP, apN, apP;
    SX_FLT alpha, beta, aii = 0;

    // indices of C-nodes
    SX_INT *cindex = (SX_INT *) sx_calloc(row, sizeof(SX_INT));

    // Step 1. Fill in values for interpolation operator P
    for (i = 0; i < row; ++i) {

        begin_row = A->Ap[i];
        end_row = A->Ap[i + 1];

        // find diagonal entry first!!!
        for (idiag = begin_row; idiag < end_row; idiag++) {
            if (A->Aj[idiag] == i) {
                aii = A->Ax[idiag];
                break;
            }
        }

        if (vec[i] == FGPT) {   // fine grid nodes
            amN = amP = apN = apP = 0.0;
            num_pcouple = 0;

            for (j = begin_row; j < end_row; ++j) {
                if (j == idiag) continue;   // skip diagonal

                // check a point strong-coupled to i or not
                IS_STRONG = FALSE;

                for (k = P->Ap[i]; k < P->Ap[i + 1]; ++k) {
                    if (P->Aj[k] == A->Aj[j]) {
                        IS_STRONG = TRUE;
                        break;
                    }
                }

                if (A->Ax[j] > 0) {
                    apN += A->Ax[j];   // sum up positive entries
                    if (IS_STRONG) {
                        apP += A->Ax[j];
                        num_pcouple++;
                    }
                }
                else {
                    amN += A->Ax[j];   // sum up negative entries
                    if (IS_STRONG) {
                        amP += A->Ax[j];
                    }
                }
            }

            // set weight factors
            alpha = amN / amP;
            if (num_pcouple > 0) {
                beta = apN / apP;
            }
            else {
                beta = 0.0;
                aii += apN;
            }

            // keep aii inside the loop to avoid floating pt error
            for (j = P->Ap[i]; j < P->Ap[i + 1]; ++j) {
                k = P->Aj[j];

                for (l = A->Ap[i]; l < A->Ap[i + 1]; l++) {
                    if (A->Aj[l] == k) break;
                }
                if (A->Ax[l] > 0) {
                    P->Ax[j] = -beta * A->Ax[l] / aii;
                }
                else {
                    P->Ax[j] = -alpha * A->Ax[l] / aii;
                }
            }
        }                       // end if vec
        else if (vec[i] == CGPT) {      // coarse grid nodes
            P->Ax[P->Ap[i]] = 1.0;
        }
    }

    // Step 2. Generate coarse level indices and set values of P.Aj
    for (index = i = 0; i < row; ++i) {
        if (vec[i] == CGPT)
            cindex[i] = index++;
    }

    P->num_cols = index;

    for (i = 0; i < P->num_nnzs; ++i) {
        j = P->Aj[i];
        P->Aj[i] = cindex[j];
    }

    // clean up
    sx_free(cindex);

    // Step 3. Truncate the prolongation operator to reduce cost
    sx_amg_interp_trunc(P, pars);
}

/**
 * \fn static void interp_STD(SX_MAT *A, SX_IVEC *vertices, SX_MAT *P, 
 *                             SX_IMAT *S, SX_AMG_PARS *pars)
 *
 * \brief Standard interpolation
 *
 * \pars A          Pointer to SX_MAT: the coefficient matrix (index starts from 0)
 * \pars vertices   Indicator vector for the C/F splitting of the variables
 * \pars P          Interpolation matrix (input: nnz pattern, output: prolongation)
 * \pars S          Strong connection matrix
 * \pars pars      Pointer to SX_AMG_PARS: AMG parseters
 *
 */
static void interp_STD(SX_MAT * A, SX_IVEC * vertices, SX_MAT * P, SX_IMAT * S, SX_AMG_PARS * pars)
{
    const SX_INT row = A->num_rows;
    SX_INT *vec = vertices->d;

    // local variables
    SX_INT i, j, k, l, m, index;
    SX_FLT alpha, factor, alN, alP;
    SX_FLT akk, akl, aik, aki;

    // indices for coarse neighbor node for every node
    SX_INT *cindex = (SX_INT *) sx_calloc(row, sizeof(SX_INT));

    // indices from column number to index in nonzeros in i-th row
    SX_INT *rindi = (SX_INT *) sx_calloc(2 * row, sizeof(SX_INT));

    // indices from column number to index in nonzeros in k-th row
    SX_INT *rindk = (SX_INT *) sx_calloc(2 * row, sizeof(SX_INT));

    // sums of strongly connected C neighbors
    SX_FLT *csum = (SX_FLT *) sx_calloc(row, sizeof(SX_FLT));

    // sums of all neighbors except ISPT
    SX_FLT *psum = (SX_FLT *) sx_calloc(row, sizeof(SX_FLT));

    // sums of all neighbors
    SX_FLT *nsum = (SX_FLT *) sx_calloc(row, sizeof(SX_FLT));

    // diagonal entries
    SX_FLT *diag = (SX_FLT *) sx_calloc(row, sizeof(SX_FLT));

    // coefficients hat a_ij for relevant CGPT of the i-th node
    SX_FLT *Ahat = (SX_FLT *) sx_calloc(row, sizeof(SX_FLT));

    // Step 0. Prepare diagonal, Cs-sum, and N-sum
    sx_iarray_set(row, cindex, -1);
    sx_blas_array_set(row, csum, 0.0);
    sx_blas_array_set(row, nsum, 0.0);

    for (i = 0; i < row; i++) {
        // set flags for strong-connected C nodes
        for (j = S->Ap[i]; j < S->Ap[i + 1]; j++) {
            k = S->Aj[j];
            if (vec[k] == CGPT) cindex[k] = i;
        }

        for (j = A->Ap[i]; j < A->Ap[i + 1]; j++) {
            k = A->Aj[j];

            if (cindex[k] == i) csum[i] += A->Ax[j];   // strong C-couplings

            if (k == i)
                diag[i] = A->Ax[j];
            else {
                nsum[i] += A->Ax[j];
                if (vec[k] != ISPT) {
                    psum[i] += A->Ax[j];
                }
            }
        }
    }

    // Step 1. Fill in values for interpolation operator P
    for (i = 0; i < row; i++) {
        if (vec[i] == FGPT) {
            alN = psum[i];
            alP = csum[i];

            // form the reverse indices for i-th row
            for (j = A->Ap[i]; j < A->Ap[i + 1]; j++) rindi[A->Aj[j]] = j;

            // clean up Ahat for relevant nodes only
            for (j = P->Ap[i]; j < P->Ap[i + 1]; j++) Ahat[P->Aj[j]] = 0.0;

            // set values of Ahat
            Ahat[i] = diag[i];

            for (j = S->Ap[i]; j < S->Ap[i + 1]; j++) {
                k = S->Aj[j];
                aik = A->Ax[rindi[k]];

                if (vec[k] == CGPT) {
                    Ahat[k] += aik;
                }
                else if (vec[k] == FGPT) {
                    akk = diag[k];

                    // form the reverse indices for k-th row
                    for (m = A->Ap[k]; m < A->Ap[k + 1]; m++) rindk[A->Aj[m]] = m;

                    factor = aik / akk;

                    // visit the strong-connected C neighbors of k, compute
                    // Ahat in the i-th row, set aki if found
                    aki = 0.0;
                    for (m = A->Ap[k]; m < A->Ap[k + 1]; m++) {
                        if (A->Aj[m] == i) {
                            aki = A->Ax[m];
                            Ahat[i] -= factor * aki;
                        }
                    }

                    for (m = S->Ap[k]; m < S->Ap[k + 1]; m++) {
                        l = S->Aj[m];
                        akl = A->Ax[rindk[l]];
                        if (vec[l] == CGPT)
                            Ahat[l] -= factor * akl;
                    }           // end for m

                    // compute Cs-sum and N-sum for Ahat
                    alN -= factor * (nsum[k] - aki + akk);
                    alP -= factor * csum[k];

                }               // end if vec[k]
            }                   // end for j

            // How about positive entries
            if (P->Ap[i + 1] > P->Ap[i]) alpha = alN / alP;

            for (j = P->Ap[i]; j < P->Ap[i + 1]; j++) {
                k = P->Aj[j];
                P->Ax[j] = -alpha * Ahat[k] / Ahat[i];
            }
        }
        else if (vec[i] == CGPT) {
            P->Ax[P->Ap[i]] = 1.0;
        }
    }                           // end for i

    // Step 2. Generate coarse level indices and set values of P.Aj
    for (index = i = 0; i < row; ++i) {
        if (vec[i] == CGPT) cindex[i] = index++;
    }

    P->num_cols = index;

    for (i = 0; i < P->Ap[P->num_rows]; ++i) {
        j = P->Aj[i];
        P->Aj[i] = cindex[j];
    }

    // clean up
    sx_free(cindex);
    sx_free(rindi);
    sx_free(rindk);
    sx_free(nsum);

    sx_free(psum);
    sx_free(csum);
    sx_free(diag);
    sx_free(Ahat);

    // Step 3. Truncate the prolongation operator to reduce cost
    sx_amg_interp_trunc(P, pars);
}
