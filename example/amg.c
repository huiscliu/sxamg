
#include "sxamg.h"
#include "mat.h"

int main(void)
{
    SX_AMG_PARS pars;
    SX_MAT A;
    SX_VEC b, x;
    char *mat_file = "A.dat";
    SX_RTN rtn;
    SX_INT prob = 3;
    SX_INT ncx = 130, ncy = 130, ncz = 14;
    SX_INT nglobal = 0;
    
    /* create distributed matrix */
    if (prob == 1) {
        nglobal = ncx;
        A = laplacian_3pt(ncx);
        sx_printf("sx: problem size: %d.\n", nglobal);
    }
    else if (prob == 2) {
        nglobal = ncx * ncx;
        A = laplacian_5pt(ncx);
        sx_printf("sx: problem size: %d, %d x %d.\n", nglobal, ncx, ncx);
    }
    else if (prob == 3) {
        nglobal = ncx * ncy * ncz;
        A = laplacian_7pt_bound(ncx, ncy, ncz);

        sx_printf("sx: problem size: %d, %d x %d x %d.\n", nglobal, ncx, ncy, ncz);
    }
    else if (prob == 4) {
        sx_mat_read(mat_file, &A);
    }
    else {
        sx_printf("sx: wrong problem.\n");
        exit(-1);
    }

    /* pars */
    sx_amg_pars_init(&pars);
    pars.maxit = 1000;
    pars.verb = 3;
    
    /* print info */
    sx_printf("A: m = %"dFMT", n = %"dFMT", nnz = %"dFMT"\n", A.num_rows,
            A.num_cols, A.num_nnzs);
    sx_amg_pars_print(&pars);

    b = sx_vec_create(A.num_rows);
    sx_vec_set_value(&b, 1.0);

    /* solve the system */
    x = sx_vec_create(A.num_rows);
    sx_vec_set_value(&x, 1.0);
    
    rtn = sx_solver_amg(&A, &x, &b, &pars);
    
    sx_printf("AMG residual: %"fFMTg"\n", rtn.ares);
    sx_printf("AMG relative residual: %"fFMTg"\n", rtn.rres);
    sx_printf("AMG iterations: %"dFMT"\n", rtn.nits);

    sx_mat_destroy(&A);
    sx_vec_destroy(&x);
    sx_vec_destroy(&b);
    
    return 0;
}
