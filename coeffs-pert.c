#include <math.h>
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

/* calculate the ket coefficients using perturbation theory */
/* see section IIIB in one-photon PRB for details */

#define COARSENESS 4   /* defines the coarseness of the step in kc with respect to
                        * the resolution that will eventually be needed for the
                        * calculation.  This can be adjusted as a tradeoff between
                        * speed and accuracy. */

#define DEBUG 0

extern Model * model;

/* timer for optimization */
GTimer * cptimer = NULL;

typedef struct {
	ThreeVector kperp;
	ThreeVector dir;
    int pos;
} odeparams;

/* define utility functions for accessing matrix elements */

#define F(m,n) fc[(m)*NB+(n)]

#define Wx(m,n) wc[(m)*NB+(n)]
#define Wy(m,n) wc[(m)*NB+(n)+N2]
#define Wz(m,n) wc[(m)*NB+(n)+2*N2]

#define O(m)    c[2*N2+m]
#define Op(m)   f[2*N2+m]

#define C(m,n) cc[(m)*NB+(n)]

static
void spew_complex_array (complex double * arr, size_t size1, size_t size2) {
    size_t m,n;
    for (m=0;m<size1;m++) {
        for (n=0;n<size2;n++) {
            printf("%1.1e,%1.1e ",creal(arr[m*size1+n]),cimag(arr[m*size1+n]));
        }
        printf("\n");
    }
    printf("\n");
}

#if DEBUG
static
void spew_matrix_complex (gsl_matrix_complex * mat) {
    size_t m,n;
    for (m=0;m<mat->size1;m++) {
        for (n=0;n<mat->size2;n++) {
            printf("%1.2e,%1.2e ",GSL_REAL(gsl_matrix_complex_get(mat, m, n)),GSL_IMAG(gsl_matrix_complex_get(mat, m, n)));
        }
        printf("\n");
    }
    printf("\n");
}
#endif

/* k.p code starts here */

/* this is the function that is evaluated at each step for the differential equation for the velocity matrix elements */

static complex double * wc;
static complex double * dHx;
static complex double * dHy;
static complex double * dHz;

static
int cfuncall (const double kc, const double c[], double f[], void *params) // calculate all components simultaneously
{
	size_t n,m,q;
	const size_t NB=model->total_bands;
    const size_t N2=NB*NB;
    
    const odeparams pars = * (odeparams *) params;
    const ThreeVector dir = pars.dir;
    
    complex double * cc = (complex double*) c;
    complex double * fc = (complex double*) f;
    
    // calculate matrix elements from c[]

    ThreeVector k = dir;
    three_vector_scale(&k,kc);
    if (!pars.pos)
        three_vector_scale(&k,-1.0);
    three_vector_incr(&k,&pars.kperp);
    
    if(fabs(dir.x[0])>1e-5*three_vector_length(&dir)) {
      model->dHdx (dHx, &k);
    }
    if(fabs(dir.x[1])>1e-5*three_vector_length(&dir)) {
      model->dHdy (dHy, &k);
    }
    if(fabs(dir.x[2])>1e-5*three_vector_length(&dir)) {
      model->dHdz (dHz, &k);
    }

    three_vector_array_project_inplace ((double*) wc, (double *) dHx, (double *) dHy, (double *) dHz, &dir, 2*N2);
    
    complex double *W = model->buffer1;
    complex_array_zero(W,N2);
    matrix_transform (W, cc, wc, NB);
    
    memset(f, 0, sizeof(double)*(2*N2+NB));

    for (n=0;n<NB;n++) {
        Op(n) = creal(W[n+n*NB]);
    }
    for (n=0;n<NB;n++) {
        for (q=0;q<NB;q++) {
            double omegadiff=O(n)-O(q);
            if (fabs(omegadiff)>1e-11) {  // was 1e-12
                complex double aa = W[q*NB+n]/omegadiff;
                for (m=0;m<NB;m++) {
                    F(m,n)+=aa*C(m,q);
                }
            }
        }
    }
    if (!pars.pos)
        for (m=0;m<2*N2+NB;m++)
        	f[m]=-f[m];
        
	return GSL_SUCCESS;
}

static void cp_init() {
  /* allocate W */
  size_t NB=model->total_bands;
  wc = complex_double_array_calloc (NB*NB);
  dHx = complex_double_array_calloc (NB*NB);
  dHy = complex_double_array_calloc (NB*NB);
  dHz = complex_double_array_calloc (NB*NB);
}

/* this function calculates the velocity matrix elements at kperp along a direction,
 * returning a list of matrix elements as a function of kc */

/* confusing but simplifies things:  kperp is in efg basis, direction in xyz basis */
/* g is the direction of the field.  See three-vector.cpp for functions that define e and f */

GList * calc_w_phi_coeffs (const ThreeVector * kperp, const ThreeVector * direction) {
    if(!cptimer) {  /* optimization timer */
        cptimer = g_timer_new();
        g_timer_start(cptimer);
    }
    else {
        g_timer_continue(cptimer);
    }
    
    cp_init();

    g_message ("starting calc_w_phi: {%1.2e, %1.2e}", kperp->x[0], kperp->x[1]);

    if (model->total_bands == 2) /* for checking the two-band PBA model */
        return calc_w_phi_constP (kperp);
    
    double te, he, kc;
    
    int status;
    size_t qq;
    
	gboolean failed = FALSE;
	
	size_t NB = model->total_bands;
	size_t N2 = (model->total_bands)*(model->total_bands);
	
	pert_ode * po = pert_ode_init (NB);
	
	// generate basis vectors e, f, g
	
	const ThreeVector *g = direction;
	const ThreeVector *e = three_vector_e (g);
	const ThreeVector *f = three_vector_f (g, e);
	
	g_message ("e is %1.2e, %1.2e, %1.2e", e->x[0], e->x[1], e->x[2]);
	g_message ("f is %1.2e, %1.2e, %1.2e", f->x[0], f->x[1], f->x[2]);
	g_message ("g is %1.2e, %1.2e, %1.2e", g->x[0], g->x[1], g->x[2]);
	
	// get kperp in basis x, y, z
	ThreeVector *kperp2 = three_vector_new (0,0,0);
	kperp2->x[0] = kperp->x[0]*e->x[0] + kperp->x[1]*f->x[0] + kperp->x[2]*g->x[0];
	kperp2->x[1] = kperp->x[0]*e->x[1] + kperp->x[1]*f->x[1] + kperp->x[2]*g->x[1];
	kperp2->x[2] = kperp->x[0]*e->x[2] + kperp->x[1]*f->x[2] + kperp->x[2]*g->x[2];
	
    if (trajectory_intersects_bad (kperp2, direction)) {
        g_timer_stop(cptimer);
        return NULL;
    }
	
	complex double * H;
	
	// set up ODE solver
	
	po->syse.function = cfuncall;
	
    odeparams params;
    po->syse.params = (void *) &params;
	
    te = 0.0;//, kc1 = KCMAX/5;
    he = 5e-7;
    
    double * call;
    complex double *callc;
	
    H = model->H(kperp2);
    double * eigval = eigenvalues_herm (H, model->total_bands);
    gsl_matrix_complex * eigvec = eigenvectors_herm (H, model->total_bands);
    m_free(H);
    size_t m, n;
    call = double_array_calloc (2*NB*NB+NB);
    callc = (complex double*) call;
    for (m=0;m<NB;m++) {
        for (n=0;n<NB;n++) {
            callc[m*NB+n]=*(complex double *)(eigvec->data + 2*(m * eigvec->tda + n));
        }
        call[2*NB*NB+m]=eigval[m];
    }
        
    d_free(eigval);
    gsl_matrix_complex_free(eigvec);
        
    params.kperp = *three_vector_copy (kperp2);
    po->syse.params = (void *) &params;
    
	params.pos = TRUE;
	params.dir = *g;
	
	/* make a copy of the initial condition because we will need it to go negative */
	
	double * cinit = double_array_clone(call,2*N2+NB);
	
	// make list of points, add first point
    
    GList * pointlist = g_list_append (NULL, matrix_element_create_from_coeffs(kperp, callc, 0.0));
	
    // evolve in positive g direction
    
    g_message("going in g direction");
	
    te = 0.0;//, kc1 = KCMAX/5;
    he = 5e-7;
    
    qq=1;
    double dkc=5e-5;
    
	for (kc=dkc*COARSENESS;kc<=KCMAX;kc=kc+dkc*COARSENESS) {  // coarse loop for ODE solver for W - positive kc
	    while (te < kc && !failed) {
	        status = pert_ode_evolve (po, &te, kc, &he, call);
            ThreeVector *kpar = three_vector_copy(direction);
            three_vector_scale(kpar,te);
            ThreeVector *kval = three_vector_add(kperp2,kpar);
            pointlist = g_list_prepend (pointlist, matrix_element_create_from_coeffs(kval, callc, te));
            
            if (qq>100000) {  /* danger with adaptive solver is it might get stuck */
	        	failed=TRUE;  /* this detects whether the number of points is getting too high */
	    	    g_warning("calc_w_phi_coeffs:  ode solver ran away, positive direction");
	    	}
	    	qq++;
	    }
	}
	
	pointlist = g_list_reverse (pointlist); // reverse list
	
	pert_ode_reset(po);  // reset the ode solver, initial condition
	
	d_free(call);
	call = cinit;
	callc = (complex double*) call;
	
	// evolve in negative g direction
	
    params.pos = FALSE;

    qq=0;
    te=0;
    he=5e-7;
    
    for (kc=dkc*COARSENESS;kc<=KCMAX;kc=kc+dkc*COARSENESS) {
	    while (te < kc && !failed) {
    	 	status = pert_ode_evolve (po, &te, kc, &he, call);
    	 	ThreeVector *kpar = three_vector_copy(direction);
            three_vector_scale(kpar,-te);
            ThreeVector *kval = three_vector_add(kperp2,kpar);
            pointlist = g_list_prepend (pointlist, matrix_element_create_from_coeffs(kval, callc, -te));
    	 	
    	 	if (qq>100000) {
	        	failed=TRUE;
	    	    g_warning("calc_w_phi_coeffs:  ode solver ran away, negative direction");
	    	}
	    	qq++;
    	}
    }
	
	pert_ode_free(po);
	
	d_free(call);
	
    g_timer_stop(cptimer);  /* optimization timer */
	
    m_free(wc);
    if(dHx)
      m_free(dHx);
    if(dHy)
      m_free(dHy);
    if(dHz)
      m_free(dHz);

	return pointlist;
}

/* for the 14-band model, we can diagonalize the Hamiltonian to find the initial condition (at kperp) */
void set_initial_condition2 (double * wall, gsl_matrix_complex * pxproj, gsl_matrix_complex * pyproj, gsl_matrix_complex * pzproj, complex double * rxzproj, complex double * ryzproj, complex double * rzzproj) {
	int m, n;
	int NB=model->total_bands;
	int N2=NB*NB;
	complex double * wc = (complex double *) wall;
		
	for (m=0;m<NB;m++) {  // set initial conditions for ODE solvers -- all bands
        for (n=0;n<NB;n++) {
            wc[m+NB*n]=_gsl_matrix_complex_get(pxproj,m,n);
            wc[m+NB*n+N2]=_gsl_matrix_complex_get(pyproj,m,n);
            wc[m+NB*n+2*N2]=_gsl_matrix_complex_get(pzproj,m,n);
            
            if (rzzproj) {
                wc[m+NB*n+3*N2]=rxzproj[m+n*NB];
                wc[m+NB*n+4*N2]=ryzproj[m+n*NB];
                wc[m+NB*n+5*N2]=rzzproj[m+n*NB];
            }
            else {
                wc[m+NB*n+3*N2]=0;
                wc[m+NB*n+4*N2]=0;
                wc[m+NB*n+5*N2]=0;
            }
        }
    }
}

/* From the D matrix, calculate the velocity matrix elements */

MatrixElement * matrix_element_create_from_coeffs (const ThreeVector *k, complex double * cc, double kcalc)
{
    size_t m;
    MatrixElement *me = matrix_element_new (k);
    me->kc = kcalc;

#if DEBUG
    printf("k={%e,%e,%e}\n", k->x[0],k->x[1],k->x[2]);
#endif
    
    size_t NB = (model->total_bands);
    size_t N2 = (model->total_bands)*(model->total_bands);
    double * c = (double*) cc;
    for (m=0;m<NB;m++) {
        me->energies[m]=c[2*N2+m];
    }
    
    // calculate matrix elements from c[]
    
    model->dHdx (dHx,k);
    model->dHdy (dHy,k);
    model->dHdz (dHz, k);
    
    matrix_transform3 (me->Wx, me->Wy, me->Wz, cc, dHx, dHy, dHz, NB);
    
    return me;
}

/* optimization timer */
double coeffs_pert_get_elapsed_time() {
    return g_timer_elapsed(cptimer, NULL);
}
