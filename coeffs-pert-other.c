#include <glib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include "pert-ode.h"
#include "matrix-element-pert.h"
#include "diagonalize.h"
#include "gsl_util.h"

extern Model * model;

/* special case for 2 band model */
/* matrix elements are constant */
/* only works for g=z */

GList * calc_w_phi_constP (ThreeVector * kperp) {
    GList * pointlist = NULL;
    
    ThreeVector * kval = three_vector_copy (kperp);
    kval->x[2]=-KCMAX;
    
    ThreeVector * kstep = three_vector_new (0,0, 5e-5);
    
    int NB = model->total_bands;
    
    double * wall = (double*) malloc((2*NB*NB*3*3+NB)*sizeof(double));  // where we will store W,R
    complex double *wallc = (complex double*) wall;
    
	int m;
	
	while (kval->x[2]<=KCMAX) {
    	complex double * H = model->H(kval);
    	double * eigval = eigenvalues_herm (H, model->total_bands);
    	
    	complex double * dhx = complex_double_array_calloc(NB*NB);
        model->dHdx (dhx, kval);
    	gsl_matrix_complex_view px = gsl_matrix_complex_view_array ((double *) dhx, NB, NB);
    	
    	complex double * dhy = complex_double_array_calloc(NB*NB);
        model->dHdy (dhy, kval);
        gsl_matrix_complex_view py = gsl_matrix_complex_view_array ((double *) dhy, NB, NB);
        
        complex double * dhz = complex_double_array_calloc(NB*NB);
        model->dHdz (dhz, kval);
        gsl_matrix_complex_view pz = gsl_matrix_complex_view_array ((double *) dhz, NB, NB);
        
        set_initial_condition2 (wall, &px.matrix, &py.matrix, &pz.matrix, NULL, NULL, NULL);
        
        m_free(dhx);
        m_free(dhy);
        m_free(dhz);
        
        MatrixElement * me = matrix_element_create(kval, wallc);
        
        me->kc = kval->x[2];
        
        pointlist = g_list_prepend (pointlist, me);
        
        /*for (m=0;m<NB*NB;m++) {
            printf("%f ",real(me->Wx[m]));
        }
        printf("\n");*/
        
        for (m=0;m<NB;m++) {
            me->energies[m]=eigval[m];
        }
        
        d_free(eigval);
        
        m_free(H);
        
        three_vector_incr (kval, kstep);
    }
    
    pointlist = g_list_reverse (pointlist); // reverse list
    //g_timer_stop(cptimer);
    
    return pointlist;
}


/* use the KKP approximation */
/* matrix elements are constant */

#if 0
GList * calc_w_phi_KKP (ThreeVector * kperp, ThreeVector * dir) {
    GList * pointlist = NULL;
    
    double * call;
    complex double *callc;
    
    /*if(!cptimer)
        cptimer = g_timer_new();
        
    g_timer_start(cptimer);*/
    
    call = handle_gammapt_degeneracy(*dir);
    callc = (complex double*) call;
    
    ThreeVector * kdir = three_vector_copy (dir);
    three_vector_scale (kdir, 1e-4);
    
    ThreeVector * kval = three_vector_copy (kperp);
    kval->x[2]=-KCMAX;
    
    ThreeVector * kstep = three_vector_new (0,0, 5e-5);
    
    size_t NB = model->total_bands;
    
    //double * wall = (double*) malloc((2*NB*NB*3*3+NB)*sizeof(double));  // where we will store W,R
    //complex double *wallc = (complex double*) wall;
    
	size_t m;
	
	while (kval->x[2]<=KCMAX) {
    	complex double * H = model->H(kval);
    	double * eigval = eigenvalues_herm (H, model->total_bands);
    	
    	MatrixElement * me = matrix_element_create_from_coeffs(kdir, callc, kval->x[2]);
    	MatrixElement * me2 = matrix_element_create_from_coeffs(kdir, callc, kval->x[2]);
        
        size_t * ind = new size_t[NB];
        gsl_sort_index (ind, me->energies, 1, NB);
        
        size_t mp, np;
        for (mp=0;mp<NB;mp++) {
    	    for (np=0;np<NB;np++) {
                size_t m=ind[mp];
                size_t n=ind[np];
                me2->Wx[mp+np*NB]=me->Wx[m+n*NB];
            }
        }
        
        /*for (m=0;m<NB*NB;m++) {
            printf("%f ",real(me->Wx[m]));
        }
        printf("\n");*/
        
        for (m=0;m<NB;m++) {
            me2->energies[m]=eigval[m];
        }
        
        pointlist = g_list_prepend (pointlist, me2);
        
        d_free(eigval);
        
        m_free(H);
        
        three_vector_incr (kval, kstep);
    }
    
    pointlist = g_list_reverse (pointlist); // reverse list
    
    /*g_timer_stop(cptimer);*/
    
    return pointlist;
}
#endif
