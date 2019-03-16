
#include "sxamg.h"

static void sx_mat_read(const char *filemat, SX_MAT * A)
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
    A->Ap = (SX_INT *) sx_mem_calloc(n + 1, sizeof(SX_INT));

    for (i = 0; i <= n; ++i) {
        fscanf(fp, "%"dFMT"\n", A->Ap + i);
    }

    nz = A->Ap[n];
    A->num_nnzs = nz;
    A->Aj = (SX_INT *) sx_mem_calloc(nz, sizeof(SX_INT));
    A->Ax = (SX_FLT *) sx_mem_calloc(nz, sizeof(SX_FLT));

    for (i = 0; i < nz; ++i) {
        fscanf(fp, "%"dFMT"\n", A->Aj + i);
    }

    for (i = 0; i < nz; ++i) fscanf(fp, "%"fFMTf"\n", A->Ax + i);

    fclose(fp);
}

int main(void)
{
    SX_AMG_PARS pars;
    SX_MAT A;
    SX_VEC b, x;
    char *datafile1 = "A.dat";
    int verb = 2;
    SX_AMG mg;
    SX_RTN rtn;
    
    /* pars */
    sx_amg_pars_init(&pars);
    pars.maxit = 1000;
    pars.verb = 2;
    
    /* system */
    sx_mat_read(datafile1, &A);
    
    /* print info */
    if (verb > 0) {
        sx_printf("\nA: m = %"dFMT", n = %"dFMT", nnz = %"dFMT"\n", A.num_rows,
                A.num_cols, A.num_nnzs);

        sx_amg_pars_print(&pars);
    }

    // Step 1: AMG setup phase
    sx_amg_setup(&mg, &A, &pars);

    // Step 2: AMG solve phase
    b = sx_vec_create(A.num_rows);
    sx_vec_set_value(&b, 1.0);
    
    x = sx_vec_create(A.num_rows);
    sx_vec_set_value(&x, 1.0);

    /* solve */
    rtn = sx_solver_amg_solve(&mg, &x, &b);

    sx_printf("AMG residual: %"fFMTg"\n", rtn.ares);
    sx_printf("AMG relative residual: %"fFMTg"\n", rtn.rres);
    sx_printf("AMG iterations: %"dFMT"\n", rtn.nits);

    sx_mat_destroy(&A);
    sx_vec_destroy(&x);
    sx_vec_destroy(&b);
    sx_amg_data_destroy(&mg);
    
    return 0;
}
