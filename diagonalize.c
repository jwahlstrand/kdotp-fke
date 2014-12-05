#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include "k-dot-p.h"
#include "diagonalize.h"
#include "gsl_util.h"

extern Model *model;

/* this code is not optimized -- it is not called very often */

/* return the sorted real eigenvalues of a hermitian matrix */
double * eigenvalues_herm ( const complex double * A , const size_t N) {
    /* the matrix is destroyed by gsl_eigen_herm(), so we copy it */
    complex double * A2 = complex_double_array_clone (A, N*N);
    gsl_matrix_complex_view Hv = gsl_matrix_complex_view_array ((double *) A2, N, N);
    
	gsl_vector * eig = gsl_vector_calloc(N);
   	gsl_eigen_herm_workspace * ws = gsl_eigen_herm_alloc (N);
    gsl_eigen_herm (&Hv.matrix, eig, ws);
    gsl_sort_vector(eig);
    gsl_eigen_herm_free(ws);
    m_free(A2);
    double * eigs = double_array_clone(gsl_vector_ptr(eig,0),N);
    gsl_vector_free(eig);
    return eigs;
}

gsl_matrix_complex * eigenvectors_herm (const complex double * A, const size_t N) {  // return the eigenvectors of H, sorted by energy
    gsl_vector * eig = gsl_vector_alloc(N);
    gsl_matrix_complex * eigvec = gsl_matrix_complex_alloc(N,N);
    gsl_eigen_hermv_workspace * ws2 = gsl_eigen_hermv_alloc (N);
    
    size_t i,j;
    gsl_matrix_complex * H = gsl_matrix_complex_alloc (N, N);
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            _gsl_matrix_complex_set(H, i, j, A[i*N+j]);
        }
    }
    
    gsl_eigen_hermv (H, eig, eigvec, ws2);
    gsl_eigen_hermv_sort (eig, eigvec, GSL_EIGEN_SORT_VAL_ASC);
    
    gsl_matrix_complex_free(H);
    gsl_vector_free(eig);
    gsl_eigen_hermv_free(ws2);
    
    return eigvec;
}

complex double * project_mat (gsl_matrix_complex * eigvec, complex double * mat) {
    size_t NB = eigvec->size1;
    gsl_matrix_complex * mat2 = mat_from_array (mat, NB, NB);
    gsl_matrix_complex * interm = gsl_matrix_complex_calloc (NB,NB);
    
    gsl_complex zero, one;
    GSL_SET_COMPLEX (&zero, 0, 0);
	GSL_SET_COMPLEX (&one,1.0,0.0);

    gsl_blas_zhemm(CblasLeft,CblasUpper, one, mat2, eigvec, zero, interm);  // matrix multiplication
    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans, one, eigvec, interm, zero, mat2);  // matrix multiplication
    
    gsl_matrix_complex_free (interm);
    complex double *nmat = array_from_mat (mat2, NB, NB);
    return nmat;
}

gsl_matrix_complex * Pproj (const gsl_matrix_complex * eigvec, const ThreeVector *k, int direction)  // project the momentum matrix using eigenvectors
{
    int NB = model->total_bands;
    gsl_matrix_complex * pproj = gsl_matrix_complex_calloc (NB,NB);
    gsl_matrix_complex * interm = gsl_matrix_complex_calloc (NB,NB);
    
    gsl_complex zero, one;
    GSL_SET_COMPLEX (&zero, 0, 0);
	GSL_SET_COMPLEX (&one,1.0,0.0);
	
	complex double * dH = complex_double_array_calloc(NB*NB);
	if (direction==0)
    	model->dHdx (dH, k);
    else if (direction==1)
        model->dHdy (dH, k);
    else
        model->dHdz (dH, k);
	
	gsl_matrix_complex * dHm = mat_from_array (dH, NB, NB);
    
    gsl_blas_zhemm(CblasLeft,CblasUpper, one, dHm, eigvec, zero, interm);  // matrix multiplication
    gsl_blas_zgemm(CblasConjTrans,CblasNoTrans, one, eigvec, interm, zero, pproj);  // matrix multiplication
    
    gsl_matrix_complex_free (interm);
    gsl_matrix_complex_free(dHm);
    
    return pproj;
}

