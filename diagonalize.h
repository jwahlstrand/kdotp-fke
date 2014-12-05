#include "config.h"
#include <complex.h>

#include <stdio.h>
#include <gsl/gsl_matrix_complex_double.h>
#include "three-vector.h"

double * eigenvalues_herm ( const complex double * A, const size_t N); // return the sorted real eigenvalues of a hermitian matrix
gsl_matrix_complex * eigenvectors_herm (const complex double * H, const size_t N);

complex double * project_mat (gsl_matrix_complex * eigvec, complex double * mat);
gsl_matrix_complex * Pproj (const gsl_matrix_complex * eigvec, const ThreeVector *k, int direction);
