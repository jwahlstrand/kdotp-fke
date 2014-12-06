#include <math.h>
#include "matrix-element-pert.h"
#include "gamma-spline.h"
#include <gsl/gsl_sort.h>

extern Model * model;
GTimer * gstimer = NULL;

/* GammaSpline */

/* take a list of matrix elements and makes interpolation objects,   *
 * which are then used to make Gamma                                 */

GammaSpline * gamma_spline_new (int size) {
    GammaSpline *gs = g_slice_new(GammaSpline);
    gs->k = double_array_calloc(size);
    gs->phi = double_array_calloc(size);
    gs->Wxr = double_array_calloc(size);
    gs->Wxi = double_array_calloc(size);
    gs->Wyr = double_array_calloc(size);
    gs->Wyi = double_array_calloc(size);
    gs->Wzr = double_array_calloc(size);
    gs->Wzi = double_array_calloc(size);
    
    gs->phi_spline = gsl_interp_alloc(gsl_interp_akima,size);
    gs->phi_accel = gsl_interp_accel_alloc();
    
    gs->wx_spline_re = gsl_interp_alloc(gsl_interp_akima,size);
    gs->w_accel = gsl_interp_accel_alloc();
    gs->wx_spline_im = gsl_interp_alloc(gsl_interp_akima,size);
    
    gs->wy_spline_re = gsl_interp_alloc(gsl_interp_akima,size);
    gs->wy_spline_im = gsl_interp_alloc(gsl_interp_akima,size);
    
    gs->wz_spline_re = gsl_interp_alloc(gsl_interp_akima,size);
    gs->wz_spline_im = gsl_interp_alloc(gsl_interp_akima,size);
    return gs;
}

/* reset the spline for different bands */

/* we want the bands in order of increasing energy */

void gamma_spline_reset (GammaSpline *gs, GList * const sd, int mp, int np, gboolean sort) {
	int q=0;
	
	int NB = model->total_bands;
	
	int m,n;
	
	GList * iter = sd;
	const int N = g_list_length(iter);
    
	MatrixElement *firstme = (MatrixElement *) iter->data;
	
	size_t ind[NB];
	
	if (sort) {
    	gsl_sort_index (ind, firstme->energies, 1, NB);
    	m=ind[mp];
    	n=ind[np];  // this breaks the 14 band model!
    }
    else {
    	m=mp;
    	n=np;
    }
		
	double *en = double_array_calloc(N);
	
    while (iter) {
        MatrixElement *me = (MatrixElement *) iter->data;
        gs->k[q]=me->kc;
        
        en[q]=(me->energies[n] - me->energies[m])*3e5/1240.7;
        
        gs->Wxr[q]=creal(me->Wx[m+n*NB]);
        gs->Wxi[q]=cimag(me->Wx[m+n*NB]);
        
        gs->Wyr[q]=creal(me->Wy[m+n*NB]);
        gs->Wyi[q]=cimag(me->Wy[m+n*NB]);
        
        gs->Wzr[q]=creal(me->Wz[m+n*NB]);
        gs->Wzi[q]=cimag(me->Wz[m+n*NB]);
        
        iter=g_list_next(iter);
        q++;
    }
	
    gsl_interp * en_spline;
    gsl_interp_accel * en_accel;
    
  	en_spline = gsl_interp_alloc(gsl_interp_akima,N);
    en_accel = gsl_interp_accel_alloc();
    
	q=gsl_interp_init(en_spline,gs->k,en,N);
    gsl_interp_accel_reset (en_accel);
    
    if (m!=n) {
        gs->phi[0]=0;
        for (q=1;q<N;q++) {
            gs->phi[q]=gs->phi[q-1]+2*M_PI*gsl_interp_eval_integ (en_spline, gs->k, en, gs->k[q-1], gs->k[q], en_accel);
        }
    }
    else {
        for (q=0;q<N;q++) {
            gs->phi[q]=0;
        }
    }
    
    d_free(en);
    
    gsl_interp_free(en_spline);
    gsl_interp_accel_free(en_accel);
	
	q=gsl_interp_init(gs->phi_spline,gs->k,gs->phi,N);
    gsl_interp_accel_reset (gs->phi_accel);
	
	q=gsl_interp_init(gs->wx_spline_re,gs->k,gs->Wxr,N);
    gsl_interp_accel_reset (gs->w_accel);
    
	q=gsl_interp_init(gs->wx_spline_im,gs->k,gs->Wxi,N);
	q=gsl_interp_init(gs->wy_spline_re,gs->k,gs->Wyr,N);
	q=gsl_interp_init(gs->wy_spline_im,gs->k,gs->Wyi,N);
	q=gsl_interp_init(gs->wz_spline_re,gs->k,gs->Wzr,N);
	q=gsl_interp_init(gs->wz_spline_im,gs->k,gs->Wzi,N);
}

double * gamma_spline_energy(GammaSpline *gs, GList * const sd, double kc)
{
    int NB = model->total_bands;
    size_t m;
    
    double * en = double_array_calloc(NB);
    
    for (m=0;m<NB;m++) {
        gamma_spline_reset(gs,sd,m,6,FALSE);
        en[m]=gsl_interp_eval_deriv(gs->phi_spline,gs->k,en,kc,gs->phi_accel)*1240.7/3e5/2/M_PI;
    }
    return en;
}

complex double * gamma_spline_V(GammaSpline *gs, GList * const sd, double kc, int component)
{
    if (!gstimer) {
        gstimer = g_timer_new();
        g_timer_start(gstimer);
    }
    else g_timer_continue(gstimer);
    
    int NB = model->total_bands;
    
    complex double * V = complex_double_array_calloc (NB*NB);
    
    size_t m, n;

    for (m=0;m<NB;m++) {
        for (n=0;n<NB;n++) {
            gamma_spline_reset(gs,sd,m,n,FALSE);
            if (component==0) {
                V[m+n*NB]=gsl_interp_eval(gs->wx_spline_re,gs->k,gs->Wxr,kc,gs->w_accel)+I*gsl_interp_eval(gs->wx_spline_im,gs->k,gs->Wxi,kc,gs->w_accel);
            }
            else if (component==1) {
                V[m+n*NB]=gsl_interp_eval(gs->wy_spline_re,gs->k,gs->Wyr,kc,gs->w_accel)+I*gsl_interp_eval(gs->wy_spline_im,gs->k,gs->Wyi,kc,gs->w_accel);
            }
            else if (component==2) {
                V[m+n*NB]=gsl_interp_eval(gs->wz_spline_re,gs->k,gs->Wzr,kc,gs->w_accel)+I*gsl_interp_eval(gs->wz_spline_im,gs->k,gs->Wzi,kc,gs->w_accel);
            }
        }
    }
    g_timer_stop(gstimer);
    return V;
}

/* use splines to create a Gamma array */
/* number of points in the array depends on epsilon */

Gamma * gamma_spline_calc_gamma (const GammaSpline *gs, double epsilon) {
	//int q=0;
	double kc;
    int Nkc = Nkc_from_epsilon(epsilon);
    Gamma * g = gamma_new (Nkc);
    double dkc=2*KCMAX/Nkc;
    double denom=pow(KCMAX,4);
    
    int q=Nkc/2;   /* we do this in a funny order to keep phase oscillations out of the spectrum */
	for (kc=-KCMAX+dkc;kc<KCMAX;kc=kc+dkc) {
	    complex double a = gsl_interp_eval(gs->wx_spline_re,gs->k,gs->Wxr,kc,gs->w_accel)+I*gsl_interp_eval(gs->wx_spline_im,gs->k,gs->Wxi,kc,gs->w_accel);
	    complex double b = cexp(-4*pow(fabs(kc),4)/denom+I*gsl_interp_eval(gs->phi_spline,gs->k,gs->phi,kc,gs->phi_accel)/epsilon);
	    g->x[q]= a*b;
	    
	    a = gsl_interp_eval(gs->wy_spline_re,gs->k,gs->Wyr,kc,gs->w_accel)+I*gsl_interp_eval(gs->wy_spline_im,gs->k,gs->Wyi,kc,gs->w_accel);
	    g->y[q]=a*b;
	    
	    a = gsl_interp_eval(gs->wz_spline_re,gs->k,gs->Wzr,kc,gs->w_accel)+I*gsl_interp_eval(gs->wz_spline_im,gs->k,gs->Wzi,kc,gs->w_accel);
	    
	    g->z[q]=a*b;
	    
	    q++;
	    if (q==Nkc)
	        q=0;
	}
	return g;
}

/* create little gamma (one-photon), which is used to calculate the zero field absorption */
/* we use N = 16384, which is an arbitrarily chosen value */
/* it seems to be enough to create a pretty smooth spectrum */

Gamma * gamma_spline_calc_little_gamma (GammaSpline *gs, double * omega) { // without the one over omega
    if (!gstimer) {
        gstimer = g_timer_new();
        g_timer_start(gstimer);
    }
    else g_timer_continue(gstimer);

    int q=0;
	double kc;
    int Nkc = 16384;
    Gamma * g = gamma_new (Nkc);
    
    double dkc=2*KCMAX/Nkc;
    double denom=pow(KCMAX,4);
	for (kc=-KCMAX+dkc;kc<KCMAX;kc=kc+dkc) {
	    complex double a = gsl_interp_eval(gs->wx_spline_re,gs->k, gs->Wxr, kc,gs->w_accel)+I*gsl_interp_eval(gs->wx_spline_im,gs->k,gs->Wxi,kc,gs->w_accel);
	    double b = gsl_interp_eval_deriv(gs->phi_spline,gs->k, gs->phi, kc,gs->phi_accel);
	    omega[q]=b;
	    double damp = exp(-4*pow(fabs(kc),4)/denom);
	    
	    g->x[q]=a*damp;
	    
	    complex double c = gsl_interp_eval(gs->wy_spline_re,gs->k,gs->Wyr,kc,gs->w_accel)+I*gsl_interp_eval(gs->wy_spline_im,gs->k,gs->Wyi,kc,gs->w_accel);
	    g->y[q]=c*damp;
	    
	    complex double d = gsl_interp_eval(gs->wz_spline_re,gs->k,gs->Wzr,kc,gs->w_accel)+I*gsl_interp_eval(gs->wz_spline_im,gs->k,gs->Wzi,kc,gs->w_accel);
	    g->z[q]=d*damp;
	    
	    q++;
	}
    g_timer_stop(gstimer);
    
	return g;
}

void gamma_spline_free(GammaSpline *gs) {
	fftw_free(gs->k);
	fftw_free(gs->phi);
	fftw_free(gs->Wxr);
	fftw_free(gs->Wxi);
	fftw_free(gs->Wyr);
	fftw_free(gs->Wyi);
	fftw_free(gs->Wzr);
	fftw_free(gs->Wzi);
	
	gsl_interp_free(gs->phi_spline);
    gsl_interp_free(gs->wx_spline_re);
    gsl_interp_free(gs->wy_spline_re);
    gsl_interp_free(gs->wz_spline_re);
    gsl_interp_free(gs->wx_spline_im);
    gsl_interp_free(gs->wy_spline_im);
    gsl_interp_free(gs->wz_spline_im);
    
    gsl_interp_accel_free(gs->phi_accel);
    gsl_interp_accel_free(gs->w_accel);
}

double gamma_spline_get_elapsed_time() {
    if (gstimer) {
        return g_timer_elapsed(gstimer, NULL);
    }
    else return 0.0;
}
