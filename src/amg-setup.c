
#include "amg-setup.h"
#include "internal.h"
#include <assert.h>

SX_AMG * sx_amg_setup(SX_MAT *A, SX_AMG_PARS *pars)
{
    SX_INT verb;
    SX_INT min_cdof;

    SX_INT m;
    SX_INT status = 0;
    SX_INT lvl = 0, max_lvls;
    SX_FLT setup_start, setup_end;
    SX_IMAT Scouple;
    SX_INT size;
    SX_AMG *mg;
    SX_IVEC vertices;

    assert(pars != NULL);
    assert(pars->max_levels > 0);
    assert(A->num_rows == A->num_cols);
    assert(A->num_rows > 0);
    assert(A->num_nnzs > 0);

    verb = pars->verb;
    min_cdof = SX_MAX(pars->coarse_dof, MIN_CDOF);
    max_lvls = pars->max_levels;

    /* timer */
    sx_gettime(&setup_start);

    /* create mg */
    mg = sx_amg_data_create(pars);

    /* init */
    m = A->num_rows;
    vertices = sx_ivec_create(m);

    /* init matrix */
    mg[0].A = sx_mat_struct_create(m, m, A->num_nnzs);
    sx_mat_cp(A, &mg[0].A);

    // Main AMG setup loop
    while ((mg[lvl].A.num_rows > min_cdof) && (lvl < max_lvls - 1)) {
        /* init */
        bzero(&Scouple, sizeof(Scouple));

        /*-- Coarsening and form the structure of interpolation --*/
        status = sx_amg_coarsen(&mg[lvl].A, &vertices, &mg[lvl].P, &Scouple, pars);

        // Check 1: Did coarsening step succeeded?
        if (status < 0) {
            /*-- Clean up Scouple generated in coarsening --*/
            sx_mem_free(Scouple.Ap);
            sx_mem_free(Scouple.Aj);

            // When error happens, stop at the current multigrid level!
            if (verb > 1) {
                sx_printf("### WARNING: Could not find any C-variables!\n");
                sx_printf("### WARNING: RS coarsening on level-%"dFMT" failed!\n", lvl);
            }
            status = 0;
            break;
        }

        // Check 2: Is coarse sparse too small?
        if (mg[lvl].P.num_cols < min_cdof) {
            /*-- Clean up Scouple generated in coarsening --*/
            sx_mem_free(Scouple.Ap);
            sx_mem_free(Scouple.Aj);

            break;
        }

        // Check 3: Does this coarsening step too aggressive?
        if (mg[lvl].P.num_rows > mg[lvl].P.num_cols * 10) {
            if (verb > 1) {
                sx_printf("### WARNING: Coarsening might be too aggressive!\n");
                sx_printf("### WARNING: Fine level = %"dFMT", coarse level = %"dFMT". Discard!\n",
                     mg[lvl].P.num_rows, mg[lvl].P.num_cols);
            }
        }

        /*-- Perform aggressive coarsening only up to the specified level --*/
        if (mg[lvl].P.num_cols * 1.5 > mg[lvl].A.num_rows) pars->cs_type = SX_COARSE_RS;

        /*-- Store the C/F marker --*/
        size = mg[lvl].A.num_rows;

        mg[lvl].cfmark = sx_ivec_create(size);
        memcpy(mg[lvl].cfmark.d, vertices.d, size * sizeof(SX_INT));

        /*-- Form interpolation --*/
        sx_amg_interp(&mg[lvl].A, &vertices, &mg[lvl].P, &Scouple, pars);

        /*-- Form coarse level matrix: two RAP routines available! --*/
        mg[lvl].R = sx_mat_trans(&mg[lvl].P);

        mg[lvl + 1].A = sx_blas_mat_rap(&mg[lvl].R, &mg[lvl].A, &mg[lvl].P);

        /*-- Clean up Scouple generated in coarsening --*/
        sx_mem_free(Scouple.Ap);
        sx_mem_free(Scouple.Aj);

        // Check 4: Is the coarse matrix too dense?
        if (mg[lvl].A.num_nnzs / mg[lvl].A.num_rows > mg[lvl].A.num_cols * 0.2) {
            if (verb > 1) {
                sx_printf("### WARNING: Coarse matrix is too dense!\n");
                sx_printf("### WARNING: m = n = %"dFMT", nnz = %"dFMT"!\n", mg[lvl].A.num_cols,
                        mg[lvl].A.num_nnzs);
            }

            /* free A */
            sx_mat_destroy(&mg[lvl + 1].A);

            break;
        }

        lvl++;
    }

    // setup total level number and current level
    mg[0].num_levels = max_lvls = lvl + 1;
    mg[0].wp = sx_vec_create(m);

    for (lvl = 1; lvl < max_lvls; ++lvl) {
        SX_INT mm = mg[lvl].A.num_rows;

        mg[lvl].num_levels = max_lvls;
        mg[lvl].b = sx_vec_create(mm);
        mg[lvl].x = sx_vec_create(mm);

        // allocate work arrays for the solve phase
        mg[lvl].wp = sx_vec_create(2 * mm);
    }

    sx_ivec_destroy(&vertices);

    if (verb > 1) {
        sx_gettime(&setup_end);
        sx_amg_complexity_print(mg);
        sx_printf("AMG setup time: %"fFMTg" s\n", setup_end - setup_start);
    }

    return mg;
}
