
#include "op-coarsen.h"
#include "internal.h"

#define LIST_HEAD -1 /**< head of the linked list */
#define LIST_TAIL -2 /**< tail of the linked list */

// Private routines for RS coarsening
static SX_INT cfsplitting_cls(SX_MAT *, SX_IMAT *, SX_IVEC *);
static SX_INT cfsplitting_clsp(SX_MAT *, SX_IMAT *, SX_IVEC *);
static SX_INT clean_ff_couplings(SX_IMAT *, SX_IVEC *, SX_INT, SX_INT);
static SX_INT compress_S(SX_IMAT *);

static void strong_couplings(SX_MAT *, SX_IMAT *, SX_AMG_PARS *);
static void form_P_pattern_dir(SX_MAT *, SX_IMAT *, SX_IVEC *, SX_INT, SX_INT);
static void form_P_pattern_std(SX_MAT *, SX_IMAT *, SX_IVEC *, SX_INT, SX_INT);

/**
 * \struct Link
 * \brief Struct for Links
 */
typedef struct
{
    SX_INT prev;
    SX_INT next;

} Link;

/**
 * \struct linked_list
 * \brief A linked list node
 *
 * \note This definition is adapted from hypre 2.0.
 */
typedef struct linked_list
{
    //! data
    SX_INT data;

    //! starting of the list
    SX_INT head;

    //! ending of the list
    SX_INT tail;

    //! next node
    struct linked_list *next_node;

    //! previous node
    struct linked_list *prev_node;

} ListElement;

typedef ListElement *LinkList; /**< linked list */

/**
 * \fn static void dispose_node( LinkList node_ptr )
 *
 * \brief Free memory space used by node_ptr
 *
 * \pars node_ptr   Pointer to the node in linked list
 *
 */
static void dispose_node(LinkList node_ptr)
{
    if (node_ptr) sx_free(node_ptr);
}

/**
 * \fn static void remove_node (LinkList *head_ptr, LinkList *tail_ptr,
 *                              SX_INT measure, SX_INT index, SX_INT *lists, SX_INT *where)
 *
 * \brief Removes a point from the lists
 *
 */
static void remove_node(LinkList * head_ptr, LinkList * tail_ptr,
        SX_INT measure, SX_INT index, SX_INT * lists, SX_INT * where)
{
    LinkList head = *head_ptr;
    LinkList tail = *tail_ptr;
    LinkList list_ptr = head;

    do {
        if (measure == list_ptr->data) {
            /* point to be removed is only point on list,
               which must be destroyed */
            if (list_ptr->head == index && list_ptr->tail == index) {
                /* removing only list, so num_left better be 0! */
                if (list_ptr == head && list_ptr == tail) {
                    head = NULL;
                    tail = NULL;
                    dispose_node(list_ptr);

                    *head_ptr = head;
                    *tail_ptr = tail;
                    return;
                }
                else if (head == list_ptr) {        /*removing 1st (max_measure) list */
                    list_ptr->next_node->prev_node = NULL;
                    head = list_ptr->next_node;
                    dispose_node(list_ptr);

                    *head_ptr = head;
                    *tail_ptr = tail;
                    return;
                }
                else if (tail == list_ptr) {        /* removing last list */
                    list_ptr->prev_node->next_node = NULL;
                    tail = list_ptr->prev_node;
                    dispose_node(list_ptr);

                    *head_ptr = head;
                    *tail_ptr = tail;
                    return;
                }
                else {
                    list_ptr->next_node->prev_node = list_ptr->prev_node;
                    list_ptr->prev_node->next_node = list_ptr->next_node;
                    dispose_node(list_ptr);

                    *head_ptr = head;
                    *tail_ptr = tail;
                    return;
                }
            }
            else if (list_ptr->head == index) { /* index is head of list */
                list_ptr->head = lists[index];
                where[lists[index]] = LIST_HEAD;
                return;
            }
            else if (list_ptr->tail == index) { /* index is tail of list */
                list_ptr->tail = where[index];
                lists[where[index]] = LIST_TAIL;
                return;
            }
            else {              /* index is in middle of list */
                lists[where[index]] = lists[index];
                where[lists[index]] = where[index];
                return;
            }
        }

        list_ptr = list_ptr->next_node;
    } while (list_ptr != NULL);

    sx_printf("### ERROR: This list is empty! %s : %d\n", __FILE__, __LINE__);
    return;
}

/**
 * \fn static LinkList create_node (SX_INT Item)
 *
 * \brief Create an node using Item for its data field
 *
 */
static LinkList create_node(SX_INT Item)
{
    LinkList new_node_ptr;

    /* Allocate memory space for the new node.
     * return with error if no space available
     */
    new_node_ptr = (LinkList) sx_calloc(1, sizeof(ListElement));
    new_node_ptr->data = Item;
    new_node_ptr->next_node = NULL;
    new_node_ptr->prev_node = NULL;
    new_node_ptr->head = LIST_TAIL;
    new_node_ptr->tail = LIST_HEAD;

    return (new_node_ptr);
}

/**
 * \fn static void enter_list (LinkList *head_ptr, LinkList *tail_ptr,
 *                             SX_INT measure, SX_INT index, SX_INT *lists, SX_INT *where)
 *
 * \brief Places point in new list
 *
 */
static void enter_list(LinkList *head_ptr, LinkList *tail_ptr,
        SX_INT measure, SX_INT index, SX_INT *lists, SX_INT *where)
{
    LinkList head = *head_ptr;
    LinkList tail = *tail_ptr;
    LinkList list_ptr;
    LinkList new_ptr;
    SX_INT old_tail;

    list_ptr = head;

    if (head == NULL) { /* no lists exist yet */
        new_ptr = create_node(measure);
        new_ptr->head = index;
        new_ptr->tail = index;
        lists[index] = LIST_TAIL;
        where[index] = LIST_HEAD;
        head = new_ptr;
        tail = new_ptr;

        *head_ptr = head;
        *tail_ptr = tail;

        return;
    }
    else {
        do {
            if (measure > list_ptr->data) {
                new_ptr = create_node(measure);

                new_ptr->head = index;
                new_ptr->tail = index;

                lists[index] = LIST_TAIL;
                where[index] = LIST_HEAD;

                if (list_ptr->prev_node != NULL) {
                    new_ptr->prev_node = list_ptr->prev_node;
                    list_ptr->prev_node->next_node = new_ptr;
                    list_ptr->prev_node = new_ptr;
                    new_ptr->next_node = list_ptr;
                }
                else {
                    new_ptr->next_node = list_ptr;
                    list_ptr->prev_node = new_ptr;
                    new_ptr->prev_node = NULL;
                    head = new_ptr;
                }

                *head_ptr = head;
                *tail_ptr = tail;

                return;
            }
            else if (measure == list_ptr->data) {
                old_tail = list_ptr->tail;
                lists[old_tail] = index;
                where[index] = old_tail;
                lists[index] = LIST_TAIL;
                list_ptr->tail = index;

                return;
            }

            list_ptr = list_ptr->next_node;
        } while (list_ptr != NULL);

        new_ptr = create_node(measure);
        new_ptr->head = index;
        new_ptr->tail = index;
        lists[index] = LIST_TAIL;
        where[index] = LIST_HEAD;
        tail->next_node = new_ptr;
        new_ptr->prev_node = tail;
        new_ptr->next_node = NULL;
        tail = new_ptr;

        *head_ptr = head;
        *tail_ptr = tail;

        return;
    }
}

/**
* \fn SX_INT sx_amg_coarsen(SX_MAT *A, SX_IVEC *vertices, SX_MAT *P,
*                           SX_IMAT *S, SX_AMG_PARS *pars)
*
* \brief Standard and aggressive coarsening schemes
*
* \pars A          Pointer to SX_MAT: Coefficient matrix (index starts from 0)
* \pars vertices   Indicator vector for the C/F splitting of the variables
* \pars P          Interpolation matrix (nonzero pattern only)
* \pars S          Strong connection matrix
* \pars pars      Pointer to SX_AMG_PARS: AMG parseters
*
* \return           0 if successed; otherwise, error information.
*
* \note vertices = 0: fine; 1: coarse; 2: isolated or special
*
*/
SX_INT sx_amg_coarsen(SX_MAT *A, SX_IVEC *vertices, SX_MAT *P, SX_IMAT *S, SX_AMG_PARS *pars)
{
    SX_COARSEN_TYPE coarse_type = pars->cs_type;
    SX_INT row = A->num_rows;
    SX_INT interp_type = pars->interp_type;
    SX_INT col = 0;

    // find strong couplings and return them in S
    strong_couplings(A, S, pars);

    switch (coarse_type) {
        case SX_COARSE_RS:        // Classical coarsening
            col = cfsplitting_cls(A, S, vertices);
            break;

        case SX_COARSE_RSP:       // Classical coarsening with positive connections
            col = cfsplitting_clsp(A, S, vertices);
            break;

        default:
            sx_exit_on_errcode(ERROR_AMG_COARSE_TYPE, __FUNCTION__);
    }

    if (col <= 0) return ERROR_UNKNOWN;

    switch (interp_type) {
        case SX_INTERP_DIR:       // Direct interpolation or ...
            col = clean_ff_couplings(S, vertices, row, col);
            form_P_pattern_dir(P, S, vertices, row, col);
            break;

        case SX_INTERP_STD:       // Standard interpolation
            form_P_pattern_std(P, S, vertices, row, col);
            break;

        default:
            sx_exit_on_errcode(ERROR_AMG_INTERP_TYPE, __FUNCTION__);
    }

    return 0;
}

/**
* \fn static void strong_couplings (SX_MAT *A, SX_IMAT *S, SX_AMG_PARS *pars)
*
* \brief Generate the set of all strong negative couplings
*
* \pars A          Coefficient matrix, the index starts from zero
* \pars S          Strong connection matrix
* \pars pars      AMG parseters
*
* \note   For flexibility, we do NOT compress S here!!! It is due to the C/F
*         splitting routines to decide when to compress S.
*
*/
static void strong_couplings(SX_MAT *A, SX_IMAT *S, SX_AMG_PARS *pars)
{
    const SX_FLT max_row_sum = pars->max_row_sum;
    const SX_FLT epsilon_str = pars->strong_threshold;
    const SX_INT row = A->num_rows, col = A->num_cols, row1 = row + 1;
    const SX_INT nnz = A->num_nnzs;

    SX_INT *ia = A->Ap, *ja = A->Aj;
    SX_FLT *aj = A->Ax;

    // local variables
    SX_INT i, j, begin_row, end_row;
    SX_FLT row_scl, row_sum;

    // get the diagonal entry of A: assume all connections are strong
    SX_VEC diag;

    diag = sx_mat_get_diag(A, 0);

    // copy the structure of A to S
    S->num_rows = row;
    S->num_cols = col;
    S->num_nnzs = nnz;
    S->Ax = NULL;
    S->Ap = (SX_INT *) sx_calloc(row1, sizeof(SX_INT));
    S->Aj = (SX_INT *) sx_calloc(nnz, sizeof(SX_INT));
    sx_iarray_cp(row1, ia, S->Ap);
    sx_iarray_cp(nnz, ja, S->Aj);

    for (i = 0; i < row; ++i) {

        // Compute row scale and row sum
        row_scl = row_sum = 0.0;
        begin_row = ia[i];
        end_row = ia[i + 1];

        for (j = begin_row; j < end_row; j++) {

            // Originally: Not consider positive entries
            // row_sum += aj[j];
            row_sum += SX_ABS(aj[j]);

            // Originally: Not consider positive entries
            // row_scl = MAX(row_scl, -aj[j]); // smallest negative
            if (ja[j] != i)
                row_scl = SX_MAX(row_scl, SX_ABS(aj[j]));   // largest abs

        }

        // Multiply by the strength threshold
        row_scl *= epsilon_str;

        // Find diagonal entries of S and remove them later
        for (j = begin_row; j < end_row; j++) {
            if (ja[j] == i) {
                S->Aj[j] = -1;
                break;
            }
        }

        // Mark entire row as weak couplings if strongly diagonal-dominant
        // Originally: Not consider positive entries
        // if ( ABS(row_sum) > max_row_sum * ABS(diag.d[i]) ) {
        if (row_sum < (2 - max_row_sum) * SX_ABS(diag.d[i])) {
            for (j = begin_row; j < end_row; j++)
                S->Aj[j] = -1;
        }
        else {
            for (j = begin_row; j < end_row; j++) {
                if (-A->Ax[j] <= row_scl) S->Aj[j] = -1;      // only n-couplings
            }
        }
    }

    sx_vec_destroy(&diag);
}

/**
* \fn static SX_INT compress_S (SX_IMAT *S)
*
* \brief Remove weak couplings from S (marked as -1)
*
* \pars S        Strong connection matrix (in: with weak, out: without weak)
*
* \return Number of cols of P
*
* \note   Compression is done in-place. Used by the C/F splitting schemes!
*/
static SX_INT compress_S(SX_IMAT * S)
{
    const SX_INT row = S->num_rows;
    SX_INT *ia = S->Ap;

    // local variables
    SX_INT index, i, j, begin_row, end_row;

    // compress S: remove weak connections and form strong coupling matrix
    for (index = i = 0; i < row; ++i) {
        begin_row = ia[i];
        end_row = ia[i + 1];

        ia[i] = index;
        for (j = begin_row; j < end_row; j++) {
            if (S->Aj[j] > -1) S->Aj[index++] = S->Aj[j];      // strong couplings
        }
    }

    S->num_nnzs = S->Ap[row] = index;

    if (S->num_nnzs <= 0) {
        return ERROR_UNKNOWN;
    }
    else {
        return 0;
    }
}

/**
* \fn static void rem_positive_ff (SX_MAT *A, SX_IMAT *Stemp, SX_IVEC *vertices)
*
* \brief Update interpolation support for positive strong couplings
*
* \pars A            Coefficient matrix, the index starts from zero
* \pars Stemp        Original strong connection matrix
* \pars vertices     Indicator vector for the C/F splitting of the variables
*
* \return Number of cols of P
*
*/
static void rem_positive_ff(SX_MAT *A, SX_IMAT *Stemp, SX_IVEC *vertices)
{
    const SX_INT row = A->num_rows;
    SX_INT *ia = A->Ap, *vec = vertices->d;

    SX_FLT row_scl, max_entry;
    SX_INT i, j, ji, max_index;

    for (i = 0; i < row; ++i) {

        if (vec[i] != FGPT)
            continue;           // skip non F-variables

        row_scl = 0.0;
        for (ji = ia[i]; ji < ia[i + 1]; ++ji) {
            j = A->Aj[ji];
            if (j == i) continue;       // skip diagonal

            row_scl = SX_MAX(row_scl, SX_ABS(A->Ax[ji]));  // max abs entry
        }                       // end for ji

        row_scl *= 0.75;

        // looking for strong F-F connections
        max_index = -1;
        max_entry = 0.0;
        for (ji = ia[i]; ji < ia[i + 1]; ++ji) {
            j = A->Aj[ji];
            if (j == i)
                continue;       // skip diagonal
            if (vec[j] != FGPT)
                continue;       // skip F-C connections

            if (A->Ax[ji] > row_scl) {
                Stemp->Aj[ji] = j;

                if (A->Ax[ji] > max_entry) {
                    max_entry = A->Ax[ji];
                    max_index = j;      // max positive entry
                }
            }
        }                       // end for ji

        // mark max positive entry as C-point
        if (max_index != -1) vec[max_index] = CGPT;
    }
}

/**
 * \fn static SX_INT cfsplitting_cls (SX_MAT *A, SX_IMAT *S, SX_IVEC *vertices)
 *
 * \brief Find coarse level variables (classic C/F splitting)
 *
 * \pars A            Coefficient matrix, the index starts from zero
 * \pars S            Strong connection matrix
 * \pars vertices     Indicator vector for the C/F splitting of the variables
 *
 * \return Number of cols of P
 *
 * \note   Coarsening Phase ONE: find coarse level points
 *
 */
static SX_INT cfsplitting_cls(SX_MAT *A, SX_IMAT *S, SX_IVEC *vertices)
{
    const SX_INT row = A->num_rows;

    // local variables
    SX_INT col = 0;
    SX_INT maxmeas, maxnode, num_left = 0;
    SX_INT measure, newmeas;
    SX_INT i, j, k, l;
    SX_INT *vec = vertices->d;
    SX_INT *work = (SX_INT *) sx_calloc(3 * row, sizeof(SX_INT));
    SX_INT *lists = work, *where = lists + row, *lambda = where + row;
    SX_IMAT ST;

    SX_INT set_empty = 1;
    SX_INT jkeep = 0, cnt, index;
    SX_INT row_end_S, ji, row_end_S_nabor, jj;
    SX_INT *graph_array = lambda;
    LinkList head = NULL, tail = NULL, list_ptr = NULL;

    // 0. Compress S and form S_transpose
    col = compress_S(S);
    if (col < 0) goto eofc;          // compression failed!!!

    ST = sx_imat_trans(S);

    // 1. Initialize lambda
    for (i = 0; i < row; ++i)
        lambda[i] = ST.Ap[i + 1] - ST.Ap[i];

    // 2. Before C/F splitting algorithm starts, filter out the variables which
    //    have no connections at all and mark them as special F-variables.
    for (i = 0; i < row; ++i) {
        if (S->Ap[i + 1] == S->Ap[i]) {
            vec[i] = ISPT;      // set i as an ISOLATED fine node
            lambda[i] = 0;
        }
        else {
            vec[i] = UNPT;      // set i as a undecided node
            num_left++;
        }
    }

    // 3. Form linked list for lambda (max to min)
    for (i = 0; i < row; ++i) {
        if (vec[i] == ISPT) continue;           // skip isolated variables

        measure = lambda[i];

        if (measure > 0) {
            enter_list(&head, &tail, lambda[i], i, lists, where);
        }
        else {
            if (measure < 0) sx_printf("### WARNING: Negative lambda[%"dFMT"]!\n", i);

            // Set variables with non-positive measure as F-variables
            vec[i] = FGPT;      // no strong connections, set i as fine node
            --num_left;

            // Update lambda and linked list after i->F
            for (k = S->Ap[i]; k < S->Ap[i + 1]; ++k) {
                j = S->Aj[k];
                if (vec[j] == ISPT)
                    continue;   // skip isolate variables
                if (j < i) {
                    newmeas = lambda[j];
                    if (newmeas > 0) {
                        remove_node(&head, &tail, newmeas, j, lists,
                                    where);
                    }
                    newmeas = ++(lambda[j]);
                    enter_list(&head, &tail, newmeas, j, lists, where);
                }
                else {
                    newmeas = ++(lambda[j]);
                }
            }
        }                       // end if measure
    }                           // end for i

    // 4. Main loop
    while (num_left > 0) {
        // pick $i\in U$ with $\max\lambda_i: C:=C\cup\{i\}, U:=U\\{i\}$
        maxnode = head->head;
        maxmeas = lambda[maxnode];

        if (maxmeas == 0) sx_printf("### WARNING: Head of the list has measure 0!\n");

        vec[maxnode] = CGPT;    // set maxnode as coarse node
        lambda[maxnode] = 0;
        --num_left;
        remove_node(&head, &tail, maxmeas, maxnode, lists, where);
        col++;

        // for all $j\in S_i^T\cap U: F:=F\cup\{j\}, U:=U\backslash\{j\}$
        for (i = ST.Ap[maxnode]; i < ST.Ap[maxnode + 1]; ++i) {
            j = ST.Aj[i];

            if (vec[j] != UNPT) continue;       // skip decided variables

            vec[j] = FGPT;      // set j as fine node
            remove_node(&head, &tail, lambda[j], j, lists, where);
            --num_left;

            // Update lambda and linked list after j->F
            for (l = S->Ap[j]; l < S->Ap[j + 1]; l++) {
                k = S->Aj[l];
                if (vec[k] == UNPT) {   // k is unknown
                    remove_node(&head, &tail, lambda[k], k, lists, where);

                    newmeas = ++(lambda[k]);
                    enter_list(&head, &tail, newmeas, k, lists, where);
                }
            }
        }                       // end for i

        // Update lambda and linked list after maxnode->C
        for (i = S->Ap[maxnode]; i < S->Ap[maxnode + 1]; ++i) {

            j = S->Aj[i];

            if (vec[j] != UNPT) continue;       // skip decided variables

            measure = lambda[j];
            remove_node(&head, &tail, measure, j, lists, where);
            lambda[j] = --measure;

            if (measure > 0) {
                enter_list(&head, &tail, measure, j, lists, where);
            }
            else {              // j is the only point left, set as fine variable
                vec[j] = FGPT;
                --num_left;

                // Update lambda and linked list after j->F
                for (l = S->Ap[j]; l < S->Ap[j + 1]; l++) {
                    k = S->Aj[l];
                    if (vec[k] == UNPT) {       // k is unknown
                        remove_node(&head, &tail, lambda[k], k, lists, where);
                        newmeas = ++(lambda[k]);
                        enter_list(&head, &tail, newmeas, k, lists, where);
                    }
                }               // end for l
            }                   // end if
        }                       // end for
    }                           // end while

    // C/F splitting of RS coarsening check C1 Criterion
    sx_iarray_set(row, graph_array, -1);
    for (i = 0; i < row; i++) {
        if (vec[i] == FGPT) {
            row_end_S = S->Ap[i + 1];
            for (ji = S->Ap[i]; ji < row_end_S; ji++) {
                j = S->Aj[ji];
                if (vec[j] == CGPT) {
                    graph_array[j] = i;
                }
            }
            cnt = 0;
            for (ji = S->Ap[i]; ji < row_end_S; ji++) {
                j = S->Aj[ji];
                if (vec[j] == FGPT) {
                    set_empty = 1;
                    row_end_S_nabor = S->Ap[j + 1];
                    for (jj = S->Ap[j]; jj < row_end_S_nabor; jj++) {
                        index = S->Aj[jj];
                        if (graph_array[index] == i) {
                            set_empty = 0;
                            break;
                        }
                    }
                    if (set_empty) {
                        if (cnt == 0) {
                            vec[j] = CGPT;
                            col++;
                            graph_array[j] = i;
                            jkeep = j;
                            cnt = 1;
                        }
                        else {
                            vec[i] = CGPT;
                            vec[jkeep] = FGPT;
                            break;
                        }
                    }
                }
            }
        }
    }

    sx_imat_destroy(&ST);

    if (head) {
        list_ptr = head;
        head->prev_node = NULL;
        head->next_node = NULL;
        head = list_ptr->next_node;
        sx_free(list_ptr);
    }

eofc:
    sx_free(work);

    return col;
}

/**
* \fn static SX_INT cfsplitting_clsp (SX_MAT *A, SX_IMAT *S, SX_IVEC *vertices)
*
 * \brief Find coarse level variables (C/F splitting with positive connections)
*
* \pars A            Coefficient matrix, the index starts from zero
* \pars S            Strong connection matrix
* \pars vertices     Indicator vector for the C/F splitting of the variables
*
* \return Number of cols of P
*
* \note   Compared with cfsplitting_cls, cfsplitting_clsp has an extra step for
*         checking strong positive couplings and pick some of them as C.
*
*/
static SX_INT cfsplitting_clsp(SX_MAT *A, SX_IMAT *S, SX_IVEC *vertices)
{
    const SX_INT row = A->num_rows;

    // local variables
    SX_INT col = 0;
    SX_INT maxmeas, maxnode, num_left = 0;
    SX_INT measure, newmeas;
    SX_INT i, j, k, l;

    SX_INT *ia = A->Ap, *vec = vertices->d;
    SX_INT *work = (SX_INT *) sx_calloc(3 * row, sizeof(SX_INT));
    SX_INT *lists = work, *where = lists + row, *lambda = where + row;

    LinkList head = NULL, tail = NULL, list_ptr = NULL;
    SX_IMAT ST;

    // 0. Compress S and form S_transpose
    SX_IMAT Stemp;
    Stemp.num_rows = S->num_rows;
    Stemp.num_cols = S->num_cols;
    Stemp.num_nnzs = S->num_nnzs;
    Stemp.Ap = (SX_INT *) sx_calloc(S->num_rows + 1, sizeof(SX_INT));
    Stemp.Aj = (SX_INT *) sx_calloc(S->num_nnzs, sizeof(SX_INT));
    sx_iarray_cp(S->num_rows + 1, S->Ap, Stemp.Ap);
    sx_iarray_cp(S->num_nnzs, S->Aj, Stemp.Aj);

    if (compress_S(S) < 0) goto eofc;          // compression failed!!!

    ST = sx_imat_trans(S);

    // 1. Initialize lambda
    for (i = 0; i < row; ++i)
        lambda[i] = ST.Ap[i + 1] - ST.Ap[i];

    // 2. Before C/F splitting algorithm starts, filter out the variables which
    //    have no connections at all and mark them as special F-variables.
    for (i = 0; i < row; ++i) {
        if ((ia[i + 1] - ia[i]) <= 1) {
            vec[i] = ISPT;      // set i as an ISOLATED fine node
            lambda[i] = 0;
        }
        else {
            vec[i] = UNPT;      // set i as a undecided node
            num_left++;
        }
    }                           // end for i

    // 3. Form linked list for lambda (max to min)
    for (i = 0; i < row; ++i) {

        if (vec[i] == ISPT)
            continue;           // skip isolated variables

        measure = lambda[i];

        if (measure > 0) {
            enter_list(&head, &tail, lambda[i], i, lists, where);
        }
        else {
            if (measure < 0)
                sx_printf("### WARNING: Negative lambda[%"dFMT"]!\n", i);

            // Set variables with non-positive measure as F-variables
            vec[i] = FGPT;      // no strong connections, set i as fine node
            --num_left;

            // Update lambda and linked list after i->F
            for (k = S->Ap[i]; k < S->Ap[i + 1]; ++k) {

                j = S->Aj[k];
                if (vec[j] == ISPT)
                    continue;   // skip isolate variables

                if (j < i) {    // only look at the previous points!!
                    newmeas = lambda[j];
                    if (newmeas > 0) {
                        remove_node(&head, &tail, newmeas, j, lists,
                                    where);
                    }
                    newmeas = ++(lambda[j]);
                    enter_list(&head, &tail, newmeas, j, lists, where);
                }
                else {          // will be checked later on
                    newmeas = ++(lambda[j]);
                }
            }
        }
    }

    // 4. Main loop
    while (num_left > 0) {
        // pick $i\in U$ with $\max\lambda_i: C:=C\cup\{i\}, U:=U\\{i\}$
        maxnode = head->head;
        maxmeas = lambda[maxnode];

        if (maxmeas == 0)
            sx_printf("### WARNING: Head of the list has measure 0!\n");

        vec[maxnode] = CGPT;    // set maxnode as coarse node
        lambda[maxnode] = 0;
        --num_left;
        remove_node(&head, &tail, maxmeas, maxnode, lists, where);
        col++;

        // for all $j\in S_i^T\cap U: F:=F\cup\{j\}, U:=U\backslash\{j\}$
        for (i = ST.Ap[maxnode]; i < ST.Ap[maxnode + 1]; ++i) {
            j = ST.Aj[i];

            if (vec[j] != UNPT)
                continue;       // skip decided variables

            vec[j] = FGPT;      // set j as fine node
            remove_node(&head, &tail, lambda[j], j, lists, where);
            --num_left;

            // Update lambda and linked list after j->F
            for (l = S->Ap[j]; l < S->Ap[j + 1]; l++) {
                k = S->Aj[l];
                if (vec[k] == UNPT) {   // k is unknown
                    remove_node(&head, &tail, lambda[k], k, lists,
                                where);
                    newmeas = ++(lambda[k]);
                    enter_list(&head, &tail, newmeas, k, lists, where);
                }
            }
        }

        // Update lambda and linked list after maxnode->C
        for (i = S->Ap[maxnode]; i < S->Ap[maxnode + 1]; ++i) {

            j = S->Aj[i];

            if (vec[j] != UNPT)
                continue;       // skip decided variables

            measure = lambda[j];
            remove_node(&head, &tail, measure, j, lists, where);
            lambda[j] = --measure;

            if (measure > 0) {
                enter_list(&head, &tail, measure, j, lists, where);
            }
            else {              // j is the only point left, set as fine variable
                vec[j] = FGPT;
                --num_left;

                // Update lambda and linked list after j->F
                for (l = S->Ap[j]; l < S->Ap[j + 1]; l++) {
                    k = S->Aj[l];
                    if (vec[k] == UNPT) {       // k is unknown
                        remove_node(&head, &tail, lambda[k], k, lists, where);
                        newmeas = ++(lambda[k]);

                        enter_list(&head, &tail, newmeas, k, lists, where);
                    }
                }
            }
        }
    }

    sx_imat_destroy(&ST);

    if (head) {
        list_ptr = head;
        head->prev_node = NULL;
        head->next_node = NULL;
        head = list_ptr->next_node;
        sx_free(list_ptr);
    }

    // Enforce F-C connections. Adding this step helps for the ExxonMobil test
    // problems! Need more tests though
    // col = clean_ff_couplings(S, vertices, row, col);

    rem_positive_ff(A, &Stemp, vertices);

    if (compress_S(&Stemp) < 0) goto eofc;          // compression failed!!!

    S->num_rows = Stemp.num_rows;
    S->num_cols = Stemp.num_cols;
    S->num_nnzs = Stemp.num_nnzs;

    sx_free(S->Ap);
    S->Ap = Stemp.Ap;
    sx_free(S->Aj);
    S->Aj = Stemp.Aj;

eofc:
    sx_free(work);

    return col;
}

/**
 * \fn static SX_INT clean_ff_couplings (SX_IMAT *S, SX_IVEC *vertices,
 *                                    SX_INT row, SX_INT col)
 *
 * \brief Clear some of the FF connections
 *
 * \pars S            Strong connection matrix
 * \pars vertices     Indicator vector for the C/F splitting of the variables
 * \pars row          Number of rows of P
 * \pars col          Number of columns of P
 *
 * \return Number of cols of P
 *
 * \note   Coarsening Phase TWO: remove some F-F connections by F->C. Need to be
 *         applied in direct and energy-min interpolations to make sure C-neighbors
 *         exist for each F-point!
 *
 */
static SX_INT clean_ff_couplings(SX_IMAT *S, SX_IVEC *vertices, SX_INT row, SX_INT col)
{
    // local variables
    SX_INT *vec = vertices->d;
    SX_INT *cindex = (SX_INT *) sx_calloc(row, sizeof(SX_INT));
    SX_INT set_empty = TRUE, C_i_nonempty = FALSE;
    SX_INT ci_tilde = -1, ci_tilde_mark = -1;
    SX_INT ji, jj, i, j, index;

    sx_iarray_set(row, cindex, -1);

    for (i = 0; i < row; ++i) {

        if (vec[i] != FGPT)
            continue;           // skip non F-variables

        for (ji = S->Ap[i]; ji < S->Ap[i + 1]; ++ji) {
            j = S->Aj[ji];
            if (vec[j] == CGPT)
                cindex[j] = i;  // mark C-neighbors
            else
                cindex[j] = -1; // reset cindex
        }

        if (ci_tilde_mark != i)
            ci_tilde = -1;      //???

        for (ji = S->Ap[i]; ji < S->Ap[i + 1]; ++ji) {

            j = S->Aj[ji];

            if (vec[j] != FGPT)
                continue;       // skip non F-variables

            // check whether there is a C-connection
            set_empty = TRUE;
            for (jj = S->Ap[j]; jj < S->Ap[j + 1]; ++jj) {
                index = S->Aj[jj];
                if (cindex[index] == i) {
                    set_empty = FALSE;
                    break;
                }
            }                   // end for jj

            // change the point i (if only F-F exists) to C
            if (set_empty) {
                if (C_i_nonempty) {
                    vec[i] = CGPT;
                    col++;
                    if (ci_tilde > -1) {
                        vec[ci_tilde] = FGPT;
                        col--;
                        ci_tilde = -1;
                    }
                    C_i_nonempty = FALSE;
                    break;
                }
                else {          // temporary set j->C and roll back
                    vec[j] = CGPT;
                    col++;
                    ci_tilde = j;
                    ci_tilde_mark = i;
                    C_i_nonempty = TRUE;
                    i--;        // roll back to check i-point again
                    break;
                }               // end if C_i_nonempty
            }                   // end if set_empty
        }                       // end for ji
    }                           // end for i

    sx_free(cindex);

    return col;
}

/**
 * \fn static void form_P_pattern_dir (SX_MAT *P, SX_IMAT *S, SX_IVEC *vertices,
 *                                     SX_INT row, SX_INT col)
 *
 * \brief Generate sparsity pattern of prolongation for direct interpolation
 *
 * \pars P         Pointer to the prolongation matrix
 * \pars S         Pointer to the set of all strong couplings matrix
 * \pars vertices  Pointer to the type of variables
 * \pars row       Number of rows of P
 * \pars col       Number of cols of P
 *
 */
static void form_P_pattern_dir(SX_MAT *P, SX_IMAT *S, SX_IVEC *vertices, SX_INT row, SX_INT col)
{
    // local variables
    SX_INT i, j, k, index;
    SX_INT *vec = vertices->d;

    // Initialize P matrix
    P->num_rows = row;
    P->num_cols = col;
    P->Ap = (SX_INT *) sx_calloc(row + 1, sizeof(SX_INT));

    // step 1: Find the structure IA of P first: using P as a counter
    for (i = 0; i < row; ++i) {
        switch (vec[i]) {
            case FGPT:         // fine grid points
                for (j = S->Ap[i]; j < S->Ap[i + 1]; j++) {
                    k = S->Aj[j];
                    if (vec[k] == CGPT)
                        P->Ap[i + 1]++;
                }
                break;

            case CGPT:         // coarse grid points
                P->Ap[i + 1] = 1;
                break;

            default:           // treat everything else as isolated
                P->Ap[i + 1] = 0;
                break;
        }
    }                           // end for i

    // Form P->Ap from the counter P
    for (i = 0; i < P->num_rows; ++i) P->Ap[i + 1] += P->Ap[i];

    P->num_nnzs = P->Ap[P->num_rows] - P->Ap[0];

    // step 2: Find the structure JA of P
    P->Aj = (SX_INT *) sx_calloc(P->num_nnzs, sizeof(SX_INT));
    P->Ax = (SX_FLT *) sx_calloc(P->num_nnzs, sizeof(SX_FLT));

    for (index = i = 0; i < row; ++i) {
        if (vec[i] == FGPT) {   // fine grid point
            for (j = S->Ap[i]; j < S->Ap[i + 1]; j++) {
                k = S->Aj[j];
                if (vec[k] == CGPT)
                    P->Aj[index++] = k;
            }                   // end for j
        }                       // end if
        else if (vec[i] == CGPT) {      // coarse grid point -- one entry only
            P->Aj[index++] = i;
        }
    }
}

/**
* \fn static void form_P_pattern_std (SX_MAT *P, SX_IMAT *S, SX_IVEC *vertices,
*                                     SX_INT row, SX_INT col)
*
* \brief Generate sparsity pattern of prolongation for standard interpolation
*
* \pars P         Pointer to the prolongation matrix
* \pars S         Pointer to the set of all strong couplings matrix
* \pars vertices  Pointer to the type of variables
* \pars row       Number of rows of P
* \pars col       Number of cols of P
*
*/
static void form_P_pattern_std(SX_MAT *P, SX_IMAT *S, SX_IVEC *vertices, SX_INT row, SX_INT col)
{
    // local variables
    SX_INT i, j, k, l, h, index;
    SX_INT *vec = vertices->d;

    // number of times a C-point is visited
    SX_INT *visited = (SX_INT *) sx_calloc(row, sizeof(SX_INT));

    P->num_rows = row;
    P->num_cols = col;
    P->Ap = (SX_INT *) sx_calloc(row + 1, sizeof(SX_INT));

    sx_iarray_set(row, visited, -1);

    // Step 1: Find the structure IA of P first: use P as a counter
    for (i = 0; i < row; ++i) {

        if (vec[i] == FGPT) {   // if node i is a F point
            for (j = S->Ap[i]; j < S->Ap[i + 1]; j++) {
                k = S->Aj[j];

                // if neighbor of i is a C point, good
                if ((vec[k] == CGPT) && (visited[k] != i)) {
                    visited[k] = i;
                    P->Ap[i + 1]++;
                }

                // if k is a F point and k is not i, look for indirect C neighbors
                else if ((vec[k] == FGPT) && (k != i)) {
                    for (l = S->Ap[k]; l < S->Ap[k + 1]; l++) { // neighbors of k
                        h = S->Aj[l];
                        if ((vec[h] == CGPT) && (visited[h] != i)) {
                            visited[h] = i;
                            P->Ap[i + 1]++;
                        }
                    }           // end for(l=S->Ap[k];l<S->Ap[k+1];l++)
                }               // end if (vec[k]==CGPT)
            }                   // end for (j=S->Ap[i];j<S->Ap[i+1];j++)
        }
        else if (vec[i] == CGPT) {      // if node i is a C point
            P->Ap[i + 1] = 1;
        }
        else {                  // treat everything else as isolated points
            P->Ap[i + 1] = 0;
        }                       // end if (vec[i]==FGPT)
    }                           // end for (i=0;i<row;++i)

    // Form P->Ap from the counter P
    for (i = 0; i < P->num_rows; ++i) P->Ap[i + 1] += P->Ap[i];

    P->num_nnzs = P->Ap[P->num_rows] - P->Ap[0];

    // Step 2: Find the structure JA of P
    P->Aj = (SX_INT *) sx_calloc(P->num_nnzs, sizeof(SX_INT));
    P->Ax = (SX_FLT *) sx_calloc(P->num_nnzs, sizeof(SX_FLT));

    sx_iarray_set(row, visited, -1);  // re-init visited array

    for (i = 0; i < row; ++i) {
        if (vec[i] == FGPT) {   // if node i is a F point
            index = 0;

            for (j = S->Ap[i]; j < S->Ap[i + 1]; j++) {
                k = S->Aj[j];

                // if neighbor k of i is a C point
                if ((vec[k] == CGPT) && (visited[k] != i)) {
                    visited[k] = i;
                    P->Aj[P->Ap[i] + index] = k;
                    index++;
                }
                // if neighbor k of i is a F point and k is not i
                else if ((vec[k] == FGPT) && (k != i)) {
                    for (l = S->Ap[k]; l < S->Ap[k + 1]; l++) { // neighbors of k
                        h = S->Aj[l];
                        if ((vec[h] == CGPT) && (visited[h] != i)) {
                            visited[h] = i;
                            P->Aj[P->Ap[i] + index] = h;
                            index++;
                        }
                    }           // end for (l=S->Ap[k];l<S->Ap[k+1];l++)
                }               // end if (vec[k]==CGPT)
            }                   // end for (j=S->Ap[i];j<S->Ap[i+1];j++)
        }
        else if (vec[i] == CGPT) {
            P->Aj[P->Ap[i]] = i;
        }
    }

    // clean up
    sx_free(visited);
}
