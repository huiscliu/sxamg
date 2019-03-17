
#ifndef SX_MMMM_H
#define SX_MMMM_H

#include "sxamg.h"
#include <assert.h>

void sx_mat_read(const char *filemat, SX_MAT * A);

SX_MAT laplacian_3pt(const SX_INT n);
SX_MAT laplacian_5pt(const SX_INT N);
SX_MAT laplacian_7pt_bound(const SX_INT Nx, const SX_INT Ny, const SX_INT Nz);

void sx_mat_read(const char *filemat, SX_MAT *A)
{
    SX_INT i, n;

    /* read the matrix from file */
    FILE *fp = fopen(filemat, "r");
    SX_INT nz;

    if (fp == NULL) {
        sx_printf("### ERROR: Cannot open %s!\n", filemat);
        exit(1);
    }

    fscanf(fp, "%"dFMT"\n", &n);
    A->num_rows = n;
    A->num_cols = n;
    A->Ap = (SX_INT *) sx_calloc(n + 1, sizeof(SX_INT));

    for (i = 0; i <= n; ++i) {
        fscanf(fp, "%"dFMT"\n", A->Ap + i);
    }

    nz = A->Ap[n];
    A->num_nnzs = nz;
    A->Aj = (SX_INT *) sx_calloc(nz, sizeof(SX_INT));
    A->Ax = (SX_FLT *) sx_calloc(nz, sizeof(SX_FLT));

    for (i = 0; i < nz; ++i) {
        fscanf(fp, "%"dFMT"\n", A->Aj + i);
    }

    for (i = 0; i < nz; ++i) fscanf(fp, "%"fFMTf"\n", A->Ax + i);

    fclose(fp);
}

SX_MAT laplacian_3pt(const SX_INT n)
{
    SX_MAT A;
    int k, left;
    int i;

    assert(n > 1);

    A.num_rows = A.num_cols = n;
    A.num_nnzs = 0;

    left = 64;
    A.Ap = sx_malloc(sizeof(SX_INT) * (n + 1));
    A.Aj = sx_malloc(sizeof(SX_INT) * left);
    A.Ax = sx_malloc(sizeof(SX_FLT) * left);

    k = 0;
    A.Ap[0] = 0;
    for (i = 0; i < A.num_rows; i++) A.Ap[i + 1] = 0;

    for (i = 0; i < n; i++) {
        if (left <= 32) {
            left += 64;

            A.Aj = sx_realloc(A.Aj, sizeof(SX_INT) * (k + left));
            A.Ax = sx_realloc(A.Ax, sizeof(SX_FLT) * (k + left));
        }

        if (i == 0) {
            A.Aj[k] = 0;
            A.Ax[k] = 2.;
            k++;
            left--;

            A.Aj[k] = 1;
            A.Ax[k] = -1;
            k++;
            left--;
        }
        else if (i == n - 1) {
            A.Aj[k] = n - 2;
            A.Ax[k] = -1;
            k++;
            left--;

            A.Aj[k] = n - 1;
            A.Ax[k] = 2;
            k++;
            left--;
        }
        else {
            A.Aj[k] = i - 1;
            A.Ax[k] = -1.;
            k++;
            left--;

            A.Aj[k] = i;
            A.Ax[k] = 2.02;
            k++;
            left--;

            A.Aj[k] = i + 1;
            A.Ax[k] = -1.;
            k++;
            left--;
        }

        A.Ap[i + 1] = k;
    }

    A.num_nnzs = k;

    return A;
}

SX_MAT laplacian_5pt(const SX_INT N)
{
    SX_MAT A;
    SX_INT nz = 0;
    SX_INT i, j;

    assert(N > 0);

    A.num_rows = N * N;
    A.num_cols = N * N;
    A.num_nnzs = 5 * N;

    A.Ap = sx_malloc(sizeof(SX_INT) * (A.num_rows + 1));
    A.Aj = sx_malloc(sizeof(SX_INT) * A.num_nnzs);
    A.Ax = sx_malloc(sizeof(SX_FLT) * A.num_nnzs);

    A.Ap[0] = 0;
    for (i = 0; i < A.num_rows; i++) A.Ap[i + 1] = 0;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            SX_INT indx = N * i + j;

            if (i > 0) {
                A.Aj[nz] = indx - N;
                A.Ax[nz] = -1;
                nz++;
            }

            if (j > 0) {
                A.Aj[nz] = indx - 1;
                A.Ax[nz] = -1;
                nz++;
            }

            A.Aj[nz] = indx;
            A.Ax[nz] = 4;
            nz++;

            if (j < N - 1) {
                A.Aj[nz] = indx + 1;
                A.Ax[nz] = -1;
                nz++;
            }

            if (i < N - 1) {
                A.Aj[nz] = indx + N;
                A.Ax[nz] = -1;
                nz++;
            }

            A.Ap[indx + 1] = nz;
        }
    }

    A.num_nnzs = nz;
    return A;
}

SX_MAT laplacian_7pt_bound(const SX_INT Nx, const SX_INT Ny, const SX_INT Nz)
{
    SX_MAT A;
    SX_INT nz = 0;
    SX_INT i, j, k;

    A.num_rows = Nx * Ny * Nz;
    A.num_cols = Nx * Ny * Nz;
    A.num_nnzs = 7 * Nx * Ny * Nz;

    A.Ap = sx_malloc(sizeof(SX_INT) * (A.num_rows + 1));
    A.Aj = sx_malloc(sizeof(SX_INT) * A.num_nnzs);
    A.Ax = sx_malloc(sizeof(SX_FLT) * A.num_nnzs);

    A.Ap[0] = 0;
    for (i = 0; i < A.num_rows; i++) A.Ap[i + 1] = 0;

    for (i = 0; i < Nz; i++) {
        for (j = 0; j < Ny; j++) {
            for (k = 0; k < Nx; k++) {
                SX_INT indx = Nx * Ny * i + Nx * j + k;

                if (i > 0) {
                    A.Aj[nz] = indx - Nx * Ny;
                    A.Ax[nz] = -1;
                    nz++;
                }

                if (j > 0) {
                    A.Aj[nz] = indx - Nx;
                    A.Ax[nz] = -1;
                    nz++;
                }

                if (k > 0) {
                    A.Aj[nz] = indx - 1;
                    A.Ax[nz] = -1;
                    nz++;
                }

                A.Aj[nz] = indx;
                A.Ax[nz] = 6;
                nz++;

                if (k < Nx - 1) {
                    A.Aj[nz] = indx + 1;
                    A.Ax[nz] = -1;
                    nz++;
                }

                if (j < Ny - 1) {
                    A.Aj[nz] = indx + Nx;
                    A.Ax[nz] = -1;
                    nz++;
                }

                if (i < Nz - 1) {
                    A.Aj[nz] = indx + Nx * Ny;
                    A.Ax[nz] = -1;
                    nz++;
                }

                A.Ap[indx + 1] = nz;
            }
        }
    }

    A.num_nnzs = nz;
    return A;
}

#endif
