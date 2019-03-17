
#include "vec.h"
#include "internal.h"
#include <assert.h>

/**
 * \fn SX_VEC sx_vec_create(SX_INT m)
 *
 * \brief Create SX_VEC data space of SX_FLT type
 *
 * \pars m    Number of rows
 *
 * \return u   The new SX_VEC
 *
 */
SX_VEC sx_vec_create(SX_INT m)
{
    SX_VEC u;

    assert(m >= 0);

    u.n = m;
    u.d = (SX_FLT *) sx_calloc(m, sizeof(SX_FLT));

    return u;
}

/**
 * \fn SX_IVEC sx_ivec_create(SX_INT m)
 *
 * \brief Create vector data space of SX_INT type
 *
 * \pars m   Number of rows
 *
 * \return u  The new SX_IVEC
 *
 */
SX_IVEC sx_ivec_create(SX_INT m)
{
    SX_IVEC u;

    assert(m >= 0);

    u.n = m;
    u.d = (SX_INT *) sx_calloc(m, sizeof(SX_INT));

    return u;
}

/**
 * \fn void sx_vec_destroy (SX_VEC *u)
 *
 * \brief Free vector data space of SX_FLT type
 *
 * \pars u   Pointer to SX_VEC which needs to be deallocated
 *
 */
void sx_vec_destroy(SX_VEC * u)
{
    if (u == NULL) return;

    sx_free(u->d);
    u->n = 0;
    u->d = NULL;
}

/**
 * \fn void sx_ivec_destroy (SX_IVEC *u)
 *
 * \brief Free vector data space of SX_INT type
 *
 * \pars u   Pointer to SX_IVEC which needs to be deallocated
 *
 * \note This function is same as sx_vec_destroy except input type.
 */
void sx_ivec_destroy(SX_IVEC * u)
{
    if (u == NULL) return;

    sx_free(u->d);
    u->n = 0;
    u->d = NULL;
}

/**
 * \fn void sx_vec_set_value(SX_VEC *x, SX_FLT val)
 *
 * \brief Initialize SX_VEC x[i]=val for i=0:n-1
 *
 * \pars n      Number of variables
 * \pars x      Pointer to SX_VEC
 * \pars val    Initial value for the vector
 *
 */
void sx_vec_set_value(SX_VEC *x, SX_FLT val)
{
    SX_INT i;
    SX_FLT *xpt = x->d;

    for (i = 0; i < x->n; ++i) xpt[i] = val;
}

/**
 * \fn void sx_ivec_set(SX_IVEC *u, SX_INT m)
 *
 * \brief Set SX_IVEC value to be m
 *
 * \pars  u    Pointer to SX_IVEC (MODIFIED)
 * \pars  m    Integer value of SX_IVEC
 *
 */
void sx_ivec_set(SX_IVEC *u, SX_INT m)
{
    SX_INT i;
    SX_INT n = u->n;

    for (i = 0; i < n; ++i) u->d[i] = m;
}

/**
 * \fn void sx_vec_cp (SX_VEC *x, SX_VEC *y) 
 *
 * \brief Copy SX_VEC x to SX_VEC y
 *
 * \pars x  Pointer to SX_VEC
 * \pars y  Pointer to SX_VEC (MODIFIED)
 *
 */
void sx_vec_cp(const SX_VEC *x, SX_VEC *y)
{
    assert(x->n > 0);

    y->n = x->n;
    memcpy(y->d, x->d, x->n * sizeof(SX_FLT));
}

SX_INT sx_vec_get_size(const SX_VEC *v)
{
    assert(v != NULL);

    return v->n;
}

void sx_vec_set_entry(SX_VEC *x, SX_INT index, SX_FLT val)
{
    assert(x != NULL);
    assert(index >= 0);
    assert(index < x->n);

    x->d[index] = val;
}

SX_FLT sx_vec_get_entry(const SX_VEC *x, SX_INT index)
{
    assert(x != NULL);
    assert(index >= 0);
    assert(index < x->n);

    return x->d[index];
}
