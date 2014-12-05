/* various useful GSL-related functions */

#include <complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_complex_math.h>
#include "k-dot-p.h"

#ifndef _INC_GSL_UTIL_H
#define _INC_GSL_UTIL_H

static inline
void _gsl_matrix_complex_hermitian_transpose (gsl_matrix_complex * mat) {
    size_t i,j;
    gsl_matrix_complex_transpose(mat);
    for (i=0;i<mat->size1;i++) {
        for (j=0;j<mat->size2;j++) {
            gsl_matrix_complex_set(mat,i,j,gsl_complex_conjugate(gsl_matrix_complex_get(mat,i,j)));
        }
    }
}

static inline
gsl_matrix_complex * mat_from_array (complex double * array, size_t size1, size_t size2) {
    size_t i,j;
    gsl_matrix_complex * m = gsl_matrix_complex_alloc (size1, size2);
    for (i=0;i<size1;i++) {
        for (j=0;j<size2;j++) {
            *(complex double *)(m->data + 2*(i * m->tda + j)) = array[i*size2+j] ;
        }
    }
    fftw_free(array);
    return m;
}

static inline
complex double * array_from_mat (gsl_matrix_complex * mat, size_t size1, size_t size2) {
    size_t i,j;
    complex double * array = complex_double_array_calloc (size1*size2);
    for (i=0;i<size1;i++) {
        for (j=0;j<size2;j++) {
            array[i*size2+j]=*(complex double *)(mat->data + 2*(i * mat->tda + j));
        }
    }
    gsl_matrix_complex_free(mat);
    return array;
}

static inline
complex double
_gsl_vector_complex_get (const gsl_vector_complex * v,
                              const size_t i)
{
  return *(complex double *)(v->data + 2*i*(v->stride));
}

static inline
void
_gsl_vector_complex_set (gsl_vector_complex * v,
                              const size_t i, complex double z)
{
  *(complex double *)(v->data + 2*i*(v->stride)) = z;
}

static inline
complex double
_gsl_matrix_complex_get (const gsl_matrix_complex * m, const size_t i, const size_t j)
{
    return *(complex double *)(m->data + 2*(i * m->tda + j)) ;
}

static inline
void
_gsl_matrix_complex_set(gsl_matrix_complex * m, 
                     const size_t i, const size_t j, const complex double x)
{
  *(complex double *)(m->data + 2*(i * m->tda + j)) = x ;
}

static inline
void
_gsl_matrix_complex_incr(gsl_matrix_complex * m, 
                     const size_t i, const size_t j, const complex double x)
{
  *(complex double *)(m->data + 2*(i * m->tda + j)) += x ;
}

void smooth (double * array, int len);

#endif

