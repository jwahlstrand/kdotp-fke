#include "gsl_util.h"
#include "pert-ode.h"
#include "diagonalize.h"
#include "k-dot-p.h"

pert_ode * pert_ode_init (size_t Nv) {
    pert_ode * po = (pert_ode *) malloc (sizeof(pert_ode));
    po->Te = gsl_odeiv_step_rkf45;
    po->se = gsl_odeiv_step_alloc (po->Te, 2*Nv*Nv+Nv);
    po->ce = gsl_odeiv_control_y_new (5e-7, 0);
    po->ee = gsl_odeiv_evolve_alloc (2*Nv*Nv+Nv);
	po->syse.jacobian = NULL;
	po->syse.dimension = 2*Nv*Nv+Nv;
	return po;
}

int pert_ode_evolve (pert_ode * po, double * te, double dk, double * he, double data[]) {
    return gsl_odeiv_evolve_apply (po->ee, po->ce, po->se, &(po->syse), te, dk, he, data);
}

void pert_ode_reset (pert_ode * po) {
    gsl_odeiv_step_reset (po->se);
	gsl_odeiv_evolve_reset (po->ee);
}

void pert_ode_free (pert_ode * po) {
    gsl_odeiv_evolve_free (po->ee);
    gsl_odeiv_control_free (po->ce);
    gsl_odeiv_step_free (po->se);
    free(po);
}

double * pert_initial_condition (complex double * mat, size_t NB) {
    size_t m, n;
    double * call = double_array_calloc (2*NB*NB+NB);  // where we will store C, omega
    
    double * eigval = eigenvalues_herm (mat, NB);
    gsl_matrix_complex * eigvec = eigenvectors_herm (mat, NB);
    
    complex double *callc = (complex double*) call;
    
    for (m=0;m<NB;m++) {
        call[2*NB*NB+m]=eigval[m];
    }
    
    for (m=0;m<NB;m++) {
        for (n=0;n<NB;n++) {
            callc[m+n*NB]=_gsl_matrix_complex_get(eigvec, m, n);
        }
    }
    
    d_free(eigval);
    gsl_matrix_complex_free(eigvec);
    
    return call;
}
