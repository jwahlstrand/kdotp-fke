#include <math.h>
#include <gsl/gsl_sort.h>
#include "k-dot-p.h"
#include "gamma-spline.h"
#include "calc-gamma.h"
#include "matrix-element-pert.h"
#include "zener.h"

extern Model * model;
GTimer * cgtimer = NULL;

/* Gamma */

/* object that represents Gamma or Gammahat for one particular band pair */

static inline
size_t next_power_of_two (const size_t input) {
	size_t N=2;
	while (input>N) {
		N=N*2;
	}
	return 2*N; //CHANGED TO HAVE ENOUGH RESOLUTION TO REMOVE BELOW GAP ARTIFACT IN TWO-PHOTON SPECTRUM
}

/* this is supposed to give the number of points required to get "just good enough" energy resolution for a given epsilon */

size_t Nkc_from_epsilon (const double epsilon) {
    if (fabs(epsilon)<1e-6)
        return 16384;
	double Emax=10.0;
	return 8*next_power_of_two(2*KCMAX*3e5*Emax/(1240.7*fabs(epsilon)));
}

void fft_init (size_t Nkc) {
    fftw_import_wisdom_from_filename("kdotp_fftw_wisdom");
    fftw_complex * dummy = (fftw_complex *) fftw_malloc (3*Nkc*sizeof(fftw_complex));
    int n = Nkc;
    fftw_plan p = fftw_plan_many_dft(1, &n, 3, dummy, NULL, 1, Nkc, dummy, NULL, 1, Nkc, FFTW_FORWARD, FFTW_MEASURE);
    fftw_destroy_plan(p);
    fftw_free(dummy);
}

/* this does the FFT, turning Gamma into Gammahat */


void gamma_do_fftw3_inplace (Gamma *g) {
	fftw_plan p;
	fftw_complex * in;
	
	in = (fftw_complex*) g->x;
	
	int n = g->Nkc;
	p = fftw_plan_many_dft(1, &n, 3, in, NULL, 1, g->Nkc, in, NULL, 1, g->Nkc, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}

void gamma_conjugate_all (Gamma *g) {
    size_t n;
    
    for (n=0;n<g->Nkc;n++) {
        g->x[n]=conj(g->x[n]);
        g->y[n]=conj(g->y[n]);
        g->z[n]=conj(g->z[n]);
    }
}

void gamma_normalize (Gamma *g, double epsilon) {
    size_t n;
    double factor = 4*KCMAX*3e5/1240.7/fabs(epsilon)/g->Nkc;
    
    for (n=0;n<g->Nkc;n++) {
        g->x[n]*=factor;
        g->y[n]*=factor;
        g->z[n]*=factor;
    }
}

size_t gamma_calculate_max_N (Gamma *g) {
    size_t n;
    double maxval = 0.0;
    for (n=0;n<g->Nkc;n++) {
        if (cabs(g->x[n])>maxval) maxval=cabs(g->x[n]);
        if (cabs(g->y[n])>maxval) maxval=cabs(g->y[n]);
        if (cabs(g->z[n])>maxval) maxval=cabs(g->z[n]);
    }
    // count backward, stop when anything is 1e-6 of the max value
    double thresh=5e-4*maxval;
    n=g->Nkc/2-1;
    size_t npos = g->Nkc/2-1;
    size_t nneg = g->Nkc/2;
    gboolean notfound=TRUE;
    while(notfound && n>0) {
        if (cabs(g->x[n])>thresh){
            npos=n;
            notfound=FALSE;
        }
        if (cabs(g->y[n])>thresh) {
            npos=n;
            notfound=FALSE;
        }
        if (cabs(g->z[n])>thresh) {
            npos=n;
            notfound=FALSE;
        }
        n--;
    }
    notfound=TRUE;
    n=g->Nkc/2;
    while(notfound && n<g->Nkc) {
        if (cabs(g->x[n])>thresh){
            nneg=n;
            notfound=FALSE;
        }
        if (cabs(g->y[n])>thresh) {
            nneg=n;
            notfound=FALSE;
        }
        if (cabs(g->z[n])>thresh) {
            nneg=n;
            notfound=FALSE;
        }
        n++;
    }
    g->max_N = MAX(npos,g->Nkc-nneg);
    return MAX(npos,g->Nkc-nneg);
}

Gamma * gamma_new (int Nk)
{
    Gamma *g = g_slice_new0(Gamma);
    g->Nkc = Nk;
    g->x = complex_double_array_calloc(3*g->Nkc);
    g->y = &g->x[g->Nkc];
    g->z = &g->x[2*g->Nkc];
    g->max_N = g->Nkc/2;
    return g;
}

void gamma_free(Gamma *g)
{
    fftw_free(g->x);
    g_slice_free(Gamma,g);
}

Gamma * gamma_new_size (Gamma * g, size_t newsize)  // newsize had better be a multiple of 2
{
    size_t n,q;
    Gamma *g2 = gamma_new(newsize);
    if (newsize<g->Nkc) {  // we are chopping down a spectrum
        for (n=0;n<g2->Nkc/2;n++) {
            g2->x[n]=g->x[n];
            g2->y[n]=g->y[n];
            g2->z[n]=g->z[n];
        }
        q=g->Nkc-1;
        for (n=g2->Nkc-1;n>=g2->Nkc/2;n--) {
            g2->x[n]=g->x[q];
            g2->y[n]=g->y[q];
            g2->z[n]=g->z[q];
            q--;
        }
    }                    // we are adding zeros in the time domain!
    else {
        for (n=0;n<(g->Nkc)/2;n++) {
            g2->x[n]=g->x[n];
            g2->y[n]=g->y[n];
            g2->z[n]=g->z[n];
            g2->x[g2->Nkc-1-n]=g->x[(g->Nkc)-1-n];
            g2->y[g2->Nkc-1-n]=g->y[(g->Nkc)-1-n];
            g2->z[g2->Nkc-1-n]=g->z[(g->Nkc)-1-n];
        }
    }
    return g2;
}

/* here, use zener code to produce all Gammas at once, then apply fft to each one to produce gammahats */

GammaHat ** calculate_gammahats (GList * matrix_elements, double epsilon, ThreeVector * dir, size_t max_N)
{
    size_t q;
	size_t NB = model->total_bands;
	
	const size_t Nkc = Nkc_from_epsilon (epsilon);
	
	if (!cgtimer) {
        cgtimer = g_timer_new();
        g_timer_start(cgtimer);
    }
    else g_timer_continue(cgtimer);
    
    fft_init(Nkc);
    
    GPtrArray * L = g_ptr_array_sized_new (NB);
    
    Gamma ** Gs;
	
	if (model->total_bands > 2) {
        GList * bundles = model->bundles;
        while(bundles) {
            Bundle * bundle = (Bundle *) bundles->data;
            
            complex double * Lsb;
            g_timer_stop(cgtimer);
            Lsb = zener_renormalize_bundle (matrix_elements, epsilon, dir, bundle);
            g_timer_continue(cgtimer);
            L->pdata[bundle->first]=Lsb;
            bundles = g_list_next(bundles);
        }
        
        size_t m,n,s,p,sp,pp;
        MatrixElement * me;
        double kc;
        const double dkc=2*KCMAX/Nkc;
        
        GList * iter = matrix_elements;
        
        Gs = (Gamma **) malloc(sizeof(Gamma*)*NB*NB);
	    for (q=0;q<NB*NB;q++) {
	        Gs[q]=NULL;
	    }
	    
	    const double denom=pow(KCMAX,4)/2;
        
        iter = matrix_elements;
        kc=-1e9;
        while ( fabs(kc)>1e-8 && iter ) { // find k=0 matrix_element
            me = (MatrixElement *) iter->data;
            kc=me->kc;
            
            iter = g_list_next (iter);
        }
        
        iter = g_list_previous (iter);
        GList * zeroelement = iter;
        
        size_t ind[NB];
        gsl_sort_index (ind, me->energies, 1, NB); //ind[p] where p is the global ordered index gives the natural index
        
        for (m=0;m<NB;m++) {
            const Bundle * bundle1 = find_bundle(m); // figure out which bundles apply here (one for m, one for n)
            
            
#if INTRA
            for (n=m;n<NB;n++) {  // loop over m and n            
#else
            for (n=m+1;n<NB;n++) {
#endif

#if UCBANDS
                if (m>7 && n>7) {}  // if this connects two upper band states, we don't need it
#else
                if (m>7 || n>7) {}  // if this involves upper band states at all, we don't need it
#endif
                else {
                Gamma * g = gamma_new(Nkc);
                
                Bundle * bundle2 = find_bundle(n);
                
                complex double *Lm = NULL;
                const size_t Nm=bundle1->list->len;
                
                complex double *Ln = NULL;
                size_t Nn=0;
                Lm = (complex double *) g_ptr_array_index(L, bundle1->first);
                size_t firstm = bundle1->first;
                size_t firstn = bundle2->first;

                Ln = (complex double *) g_ptr_array_index(L, bundle2->first);
                Nn = bundle2->list->len;
                
                iter = zeroelement;
                
                MatrixElement *next = (MatrixElement*) g_list_next(iter)->data;
                
                /* collect the elements we'll use into a more compact matrix here? */
                
                q=0;
                
                MatrixElement * meintrp = NULL;
                
                for (kc=0;kc<KCMAX;kc=kc+dkc) {  // calculate Gamma using L
                    if (kc>next->kc) {  // if we have gone beyond the next matrix element's kc, shift over until this is no longer true
                        while(kc>next->kc) {
                            iter = g_list_next (iter);
                            if (iter) {
                                me = next;
                                next = (MatrixElement *) iter->data;
                            }
                        }
                    }
                    
                    /* interpolation */
                    
                    double kfrac;
                    if (fabs(next->kc-me->kc)>1e-9)
                        kfrac = (kc-(me->kc))/((next->kc)-(me->kc));
                    else
                        kfrac = 1.0;
                    
                    double damp = exp(-pow(fabs(kc),4)/denom);
                    complex double LConj,litlx,litly,litlz;
                    
                    for (p=0;p<Nn;p++) {
            	        pp=ind[p+firstn];
            	        LConj = conj(Ln[q*Nn*Nn+p*Nn+(n-firstn)]);
            	        litlx = 0.0;
            	        litly = 0.0;
            	        litlz = 0.0;
                        for (s=0;s<Nm;s++) {
                            sp=ind[s+firstm]; // translate sorted number from bundle to index of the matrix element
                            litlx+=((1-kfrac)*me->Wx[sp+pp*NB]+kfrac*next->Wx[sp+pp*NB])*Lm[q*Nm*Nm+s*Nm+(m-firstm)];
                            litly+=((1-kfrac)*me->Wy[sp+pp*NB]+kfrac*next->Wy[sp+pp*NB])*Lm[q*Nm*Nm+s*Nm+(m-firstm)];
                            litlz+=((1-kfrac)*me->Wz[sp+pp*NB]+kfrac*next->Wz[sp+pp*NB])*Lm[q*Nm*Nm+s*Nm+(m-firstm)];
                        }
                        g->x[q]+=litlx*LConj;
                        g->y[q]+=litly*LConj;
                        g->z[q]+=litlz*LConj;
                    }
                    
                    g->x[q]*=damp;
                    g->y[q]*=damp;
                    g->z[q]*=damp;
                    
                    q++;
                }
                
                q=Nkc-1;
                
                iter=zeroelement;
            	iter=g_list_previous (iter);
            	me = (MatrixElement *) iter->data;
            	
            	GList * niter = g_list_previous(iter);
            	next = (MatrixElement *) niter->data;
            	
            	meintrp=NULL;
                
                for (kc=-dkc;kc>-KCMAX;kc=kc-dkc) {
                    while (kc<me->kc) {
                        iter = g_list_previous (iter);
                        next = me;
                        if (iter)
                            me = (MatrixElement *) iter->data;
                    }
                    
                    /* interpolation */
                    
                    double kfrac;
                    if (fabs(next->kc-me->kc)>1e-7)
                        kfrac = 1.0-(kc-(me->kc))/((next->kc)-(me->kc));
                    else
                        kfrac = 0.0;
                        
                    double damp = exp(-pow(fabs(kc),4)/denom);
                    
                    complex double LConj,litlx,litly,litlz;
                    
                    for (p=0;p<Nn;p++) {
            	        pp=ind[p+firstn];
            	        LConj = conj(Ln[q*Nn*Nn+p*Nn+(n-firstn)]);
            	        litlx = 0.0;
            	        litly = 0.0;
            	        litlz = 0.0;
                        for (s=0;s<Nm;s++) {
                            sp=ind[s+firstm]; // translate sorted number from bundle to index of the matrix element
                            litlx+=((1-kfrac)*me->Wx[sp+pp*NB]+kfrac*next->Wx[sp+pp*NB])*Lm[q*Nm*Nm+s*Nm+(m-firstm)];
                            litly+=((1-kfrac)*me->Wy[sp+pp*NB]+kfrac*next->Wy[sp+pp*NB])*Lm[q*Nm*Nm+s*Nm+(m-firstm)];
                            litlz+=((1-kfrac)*me->Wz[sp+pp*NB]+kfrac*next->Wz[sp+pp*NB])*Lm[q*Nm*Nm+s*Nm+(m-firstm)];
                        }
                        g->x[q]+=litlx*LConj;
                        g->y[q]+=litly*LConj;
                        g->z[q]+=litlz*LConj;
                    }
                    
                    g->x[q]*=damp;
                    g->y[q]*=damp;
                    g->z[q]*=damp;
                    
                    q--;
                }
                
                // pad with zeros for increased resolution
                
                Gamma * g2 = gamma_new_size(g,2*g->Nkc);
                gamma_free(g);
                                
                gamma_do_fftw3_inplace(g2); // take FFT
                gamma_normalize(g2,epsilon);
                
                // add to GammaHat
                
                Gs[n+NB*m]=gamma_new_size(g2,max_N);
                gamma_free(g2);
                
                }
            }
        }
        bundles = model->bundles;
        while(bundles) {
            Bundle * bundle = (Bundle *) bundles->data;
            
            m_free((complex double* )L->pdata[bundle->first]);
            bundles = g_list_next(bundles);
        }
        g_ptr_array_free(L, TRUE);
	}
	else {
	    Gs = (Gamma **) malloc(sizeof(Gamma*)*NB*NB);
	    for (q=0;q<NB*NB;q++) {
	        Gs[q]=NULL;
	    }
	}
		
	GammaSpline * gs;
	if (model->total_bands==2 ) {
    	gs = gamma_spline_new (g_list_length(matrix_elements));
	}
	
	Gamma * g;
	
	size_t v,c;
	for (v=0;v<NB;v++) {
#if INTRA
	    for (c=v;c<NB;c++) {  /* for off-diagonal elements, only store one direction, the other one is the complex conjugate */
#else
        for (c=v+1;c<NB;c++) {
#endif
			if (epsilon<0) {
			    gamma_conjugate_all(Gs[c+NB*v]);
			}
			
			if (model->total_bands > 2) {
			}
			else {
			    gamma_spline_reset (gs,matrix_elements, v, c, TRUE);
                g = gamma_spline_calc_gamma (gs,epsilon);
			    gamma_do_fftw3_inplace (g);
			    gamma_normalize(g,epsilon);
			    
			    Gs[c+NB*v]=g;
			}
		}
	}
	
	if (model->total_bands==2 ) {
    	gamma_spline_free(gs);
    }
    
    g_timer_stop(cgtimer);

	return Gs;
}

Bundle * find_bundle(size_t m) {
    GList * bundles = model->bundles;
    int n;
    while(bundles) {
        Bundle * bundle = (Bundle *) bundles->data;
        for (n=0;n<bundle->list->len;n++) {
            if (m==g_array_index(bundle->list,gint,n))
                return bundle;
        }
        bundles = g_list_next(bundles);
    }
    return NULL;
}

double calc_gamma_get_elapsed_time() {
    if (cgtimer) {
        return g_timer_elapsed(cgtimer, NULL);
    }
    else return 0.0;
}
