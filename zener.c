#define _USE_MATH_DEFINES
#include <gsl/gsl_odeiv.h>
#include <math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include "gsl_util.h"
#include "pert-ode.h"
#include "three-vector.h"
#include "diagonalize.h"
#include "matrix-element-pert.h"
#include "zener.h"

#define DEBUG 0

/* find the evolution matrix m */
/* most notation is consistent with John's ip.pdf */
/* which discusses the interaction picture and also the "other picture" */

extern Model * model;
extern Material * material;

GTimer * ztimer = NULL;

#define NB (model->total_bands)

static gsl_interp_accel *aTS;

#define F(m,n) fc[(m)*NBANDS+(n)]
#define C(m,n) cc[(m)*NBANDS+(n)]

typedef struct {
	double * kraw;
	double * TSraw;
	size_t N;
	gsl_interp ** TS;
	gsl_spline ** M;
	double epsilon;
	size_t NBANDS;
} dUdtparams;

static    /* this is the stepping function for the interaction picture */
int dUdtfunc_interaction (const double kc, const double c[], double f[], void *params) {
    // M is in params named "Vi"
    
    size_t m,n,p;
    dUdtparams pars = * (dUdtparams *) params;
    
    size_t N = pars.N;
    size_t NBANDS = pars.NBANDS;
    
    gsl_spline ** sM = pars.M;
    
    complex double * fc = (complex double *) f;
    complex double * cc = (complex double *) c;
    
    complex double eiM[NBANDS];
    complex double TS[NBANDS*NBANDS];
    
    const double oneovereps = 1.0/pars.epsilon;
    
    for (m=0;m<NBANDS;m++) {
        eiM[m]=cexp(I*gsl_spline_eval (sM[m], kc, aTS)*oneovereps);
        TS[m*NBANDS+m]=gsl_interp_eval (pars.TS[m*NBANDS+m], pars.kraw, &pars.TSraw[m*NBANDS*2*N+m*2*N], kc, aTS);
        for (n=m+1;n<NBANDS;n++) {
            TS[m*NBANDS+n]=gsl_interp_eval (pars.TS[m*NBANDS+n], pars.kraw, &pars.TSraw[m*NBANDS*2*N+n*2*N], kc, aTS)+I*gsl_interp_eval (pars.TS[m*NBANDS+n+NBANDS*NBANDS], pars.kraw, &pars.TSraw[m*NBANDS*2*N+n*2*N+N], kc, aTS);
        }
    }
    
    for (m=0;m<NBANDS;m++) {
        for (n=0;n<m;n++) {
            TS[m*NBANDS+n]=conj(TS[n*NBANDS+m]);
        }
    }
    
    /*complex double d;
    for (m=0;m<NBANDS;m++) {
        complex double iee = -I*eiM[m]*oneovereps;
        for (n=0;n<NBANDS;n++) {
            F(m,n)=0.0;
            for (p=0;p<NBANDS;p++) {
                d=TS[m*NBANDS+p];
                if (p==m) {
                    d-=creal(TS[m*NBANDS+m]);
                }
                F(m,n)+=d*conj(eiM[p])*C(p,n);
            }
            F(m,n)*=iee;
        }
    }*/
    /*for (m=0;m<NBANDS;m++) {
        complex double iee = -I*eiM[m]*oneovereps;
        for (n=0;n<NBANDS;n++) {
            F(m,n)=-creal(TS[m*NBANDS+m])*conj(eiM[m])*C(m,n);
            for (p=0;p<NBANDS;p++) {
                F(m,n)+=TS[m*NBANDS+p]*conj(eiM[p])*C(p,n);
            }
            F(m,n)*=iee;
        }
    }*/
    complex double d1, d2;
    for (m=0;m<NBANDS;m++) {
        complex double iee = -I*eiM[m]*oneovereps;
        for (n=0;n<NBANDS;n+=2) {
            d1=-creal(TS[m*NBANDS+m])*conj(eiM[m])*C(m,n);
            d2=-creal(TS[m*NBANDS+m])*conj(eiM[m])*C(m,n+1);
            for (p=0;p<NBANDS;p++) {
                complex double ce = conj(eiM[p])*TS[m*NBANDS+p];
                d1+=ce*C(p,n);
                d2+=ce*C(p,n+1);
            }
            F(m,n)=iee*d1;
            F(m,n+1)=iee*d2;
        }
    }
    /*complex double d1, d2, d3, d4;
    for (m=0;m<NBANDS;m+=2) {
        complex double iee1 = -I*eiM[m]*oneovereps;
        complex double iee2 = -I*eiM[m+1]*oneovereps;
        for (n=0;n<NBANDS;n+=2) {
            d1=-creal(TS[m*NBANDS+m])*conj(eiM[m])*C(m,n);
            d2=-creal(TS[m*NBANDS+m])*conj(eiM[m])*C(m,n+1);
            d3=-creal(TS[(m+1)*NBANDS+m+1])*conj(eiM[m+1])*C(m+1,n);
            d4=-creal(TS[(m+1)*NBANDS+m+1])*conj(eiM[m+1])*C(m+1,n+1);
            for (p=0;p<NBANDS;p++) {
                complex double ce1 = conj(eiM[p])*TS[m*NBANDS+p];
                d1+=ce1*C(p,n);
                d2+=ce1*C(p,n+1);
                complex double ce2 = conj(eiM[p])*TS[(m+1)*NBANDS+p];
                d3+=ce2*C(p,n);
                d4+=ce2*C(p,n+1);
            }
            F(m,n)=iee1*d1;
            F(m,n+1)=iee1*d2;
            F(m+1,n)=iee2*d3;
            F(m+1,n+1)=iee2*d4;
        }
    }*/
    
    return GSL_SUCCESS;
}

/* see if we are at a point of degeneracy by examining the band energies */

static gboolean degenerate_point (MatrixElement * me) {
    size_t mpp,npp;
    const double danger = 5e-6;
    if (model->total_bands>=14) {
    for (mpp=0;mpp<NB;mpp++) {
        for (npp=mpp+1;npp<NB;npp++) {
            if ((fabs(me->energies[mpp]-me->energies[npp])<danger)) {
                return TRUE;
            }
        }
    }
    }
    return FALSE;
}

/* take the matrix elements and calculate the renormalized Gamma for a given field strength */
/* returns a array containing L subblocks, one for each kc */

complex double * zener_renormalize_bundle (GList * list, double epsilon, ThreeVector * dir, Bundle * bundle)
{
    if(!ztimer) {
        ztimer = g_timer_new();
        g_timer_start(ztimer);
    }
    else {
        g_timer_continue(ztimer);
    }

    size_t m, n, q, n1, n2;
    size_t N = g_list_length(list);
    
    const size_t NSB = bundle->list->len;
    
    /* construct T and S splines from the matrix elements */
    
    gsl_interp * TS[2*NSB*NSB];
    aTS = gsl_interp_accel_alloc ();
    gsl_interp_accel_reset (aTS);
    
    for (n=0;n<2*NSB*NSB;n++) {
        TS[n] = gsl_interp_alloc (gsl_interp_akima, N);
    }
    
    double * k = (double *) malloc(sizeof(double)*N);
    
    GList * iter = list;
    
    MatrixElement *k0 = (MatrixElement *) iter->data;
    size_t ind[NB];
    gsl_sort_index (ind, k0->energies, 1, NB); // ind[m], where m is the global ordered index, gives the natural index
    
    double * TSraw = double_array_calloc (2*NSB*NSB*N);
    
    q=0;
    while (iter) {   /* go through list of W matrices, construct arrays that will be used to interpolate the matrix elements */
        MatrixElement *me = (MatrixElement *) iter->data;
        k[q]=me->kc;
        
        for (n1=0;n1<NSB;n1++) {
            m=ind[g_array_index(bundle->list,gint,n1)];  // need the natural index to find correct matrix element for a given band
            for (n2=0;n2<NSB;n2++) {
                n=ind[g_array_index(bundle->list,gint,n2)];
                if (n1==n2) { /* this is the diagonal part, T */
                    TSraw[n1*NSB*2*N+n2*2*N+q]=(me->energies[m])*2*M_PI*3e5/1240.7;
                    TSraw[n1*NSB*2*N+n2*2*N+N+q]=0;
                }
                else if (fabs(me->energies[m]-me->energies[n])>1e-10) {  /* the non-diagonal part, S */
                    /*  here, we need to use the direction of the DC field */
                        complex double element = -(me->Wx[m*NB+n]*dir->x[0]+me->Wy[m*NB+n]*dir->x[1]+me->Wz[m*NB+n]*dir->x[2])/(me->energies[m]-me->energies[n])/I*epsilon;
                        TSraw[n1*NSB*2*N+n2*2*N+q]=creal(element);
                        TSraw[n1*NSB*2*N+n2*2*N+N+q]=cimag(element);
                }
            }
        }
        
        iter=g_list_next(iter);
        q++;
    }
    
    for (n1=0;n1<NSB;n1++) {   /* initialize splines for interpolation */
        for (n2=n1;n2<NSB;n2++) {
            q=gsl_interp_init(TS[n1*NSB+n2],k,&TSraw[n1*NSB*2*N+n2*2*N],N);
            q=gsl_interp_init(TS[n1*NSB+n2+NSB*NSB],k,&TSraw[n1*NSB*2*N+n2*2*N+N],N);
        }
    }

    /* initialize ODE solver */
    
    const gsl_odeiv_step_type * Te;
    gsl_odeiv_step * se;
    gsl_odeiv_control * ce;
    gsl_odeiv_evolve * ee;
    gsl_odeiv_system syse;
    
    Te = gsl_odeiv_step_rkf45;
    se = gsl_odeiv_step_alloc (Te, 2*NSB*NSB);
    ce = gsl_odeiv_control_y_new (1e-6, 0);
    ee = gsl_odeiv_evolve_alloc (2*NSB*NSB);
	syse.function = dUdtfunc_interaction;
	syse.jacobian = NULL;
	syse.dimension = 2*NSB*NSB;
	
	dUdtparams params;
	params.epsilon = epsilon;
	
	/* find kpar=0 velocity matrix */
	
	iter = list;
    MatrixElement *me = NULL;
    double kc=-1e9;
    while ( fabs(kc)>1e-8 && iter ) {
        me = (MatrixElement *) iter->data;
        kc=me->kc;
        
        iter = g_list_next (iter);
    }
    
    iter = g_list_previous (iter);
    GList * zeroelement = iter;
	
	/* create splines for M */
	
	double * Mraw = double_array_calloc (N);
	gsl_spline * sM[NSB];
    
	for (m=0;m<NSB;m++) {
	    for (q=0;q<N;q++) {
            Mraw[q]=gsl_interp_eval(TS[m+NSB*m], k, &TSraw[m*NSB*2*N+m*2*N], k[q], aTS);
	    }

	    sM[m] = gsl_spline_alloc (gsl_interp_akima, N);
        q=gsl_spline_init(sM[m],k,Mraw,N);
	}
	
    /* we need the integral of M */
    
    double * iM = (double *) malloc(sizeof(double)*N);
    for (m=0;m<NSB;m++) {        
        iM[0]=0;
        for (q=1;q<N;q++) {
            iM[q]=iM[q-1]+gsl_spline_eval_integ (sM[m], k[q-1], k[q], aTS);
        }
        
        gsl_spline_free(sM[m]);
        sM[m] = gsl_spline_alloc(gsl_interp_akima,N);
        q=gsl_spline_init(sM[m],k,iM,N);
    }
    free(iM);
    
	params.M = sM;
	params.TS = TS;
	params.kraw = k;
	params.TSraw = TSraw;
	params.N = N;
	params.NBANDS = NSB;
	
	d_free(Mraw);
	
    complex double * U = complex_double_array_calloc (NSB*NSB);
	complex double * U0 = U;
	
	for (m=0;m<NSB;m++) {
	    U0[m+m*NSB]=1.0;  /* initial condition is identity matrix */
	}
	
	complex double * mU = complex_double_array_clone (U0, NSB*NSB);
	double * rU=(double *) mU;  /* real version of mU */
	
	syse.params = (void *) &params;
	
	double he, te;
	int failed, status;
	failed=FALSE;
	he=5e-8;
    te=0;
    q=1;
    
    size_t Nkc = Nkc_from_epsilon (epsilon);
    double dkc = 2*KCMAX/Nkc;
    
    complex double * Lsb = complex_double_array_calloc(Nkc*NSB*NSB); /* we output a big array contaning a bunch of L's */
    
    iter = list;
    kc=-1e9;
    while ( fabs(kc)>1e-8 && iter ) {
        me = (MatrixElement *) iter->data;
        kc=me->kc;
        
        iter = g_list_next (iter);
    }
    
    iter = g_list_previous (iter);
    
    q=0;
    
    complex double phi[NSB];
    complex_array_zero (phi, NSB);
    
    for (m=0;m<NSB;m++) {
        phi[m]=gsl_spline_eval (sM[m], 0, aTS)/epsilon;
    }
	
	complex double * tU = complex_double_array_calloc (NSB*NSB); /* we use this to make a copy of mU below */
	
	size_t SUBDIV = 1;  // substeps to take (to handle places where the matrix elements change very quickly)
    double kcp = 0.0;
    
    for (kc=0;kc<KCMAX;kc=kc+dkc) {  // coarse loop
        if (kc>0)
        for (kcp=kc-dkc;kcp<=kc;kcp=kcp+dkc/SUBDIV) {  // fine loop - not used if SUBDIV == 1
            while (te < kcp && !failed) {
                status = gsl_odeiv_evolve_apply (ee, ce, se, &syse, &te, kcp, &he, rU);  /* take a step */
            }
        }
        
        complex_double_array_copy (tU, mU, NSB*NSB);
        
        /* find phase part using M */
        
        if (kc>0) {
            for (m=0;m<NSB;m++) {
                phi[m]=gsl_spline_eval (sM[m], kc, aTS)/epsilon;
            }
        }
        
        for (m=0;m<NSB;m++) {  /* use the copy to find the actual m matrix */
            complex double phase = cexp(-I*(phi[m]));
            for (n=0;n<NSB;n++) {
                tU[m*NSB+n]*=phase;
            }
        }
        
        complex_double_array_copy (&Lsb[NSB*NSB*q],tU,NSB*NSB);  /* copy that matrix into enormous output array */
                        
        q++;
    }
    
    complex_double_array_copy (mU, U0, NSB*NSB);
    te = 0;
    he = -5e-8;
    gsl_odeiv_step_reset (se);
	gsl_odeiv_evolve_reset (ee);
	
	q=Nkc-1;
	
	iter=zeroelement;
	iter=g_list_previous (iter);
	me = (MatrixElement *) iter->data;
		
	complex_array_zero (phi, NSB);
	
    /* now go in the negative direction */
    
    for (kc=-dkc;kc>-KCMAX;kc=kc-dkc) {
        for (kcp=kc+dkc;kcp>=kc;kcp=kcp-dkc/SUBDIV) {
            while (te > kcp && !failed) {
                status = gsl_odeiv_evolve_apply (ee, ce, se, &syse, &te, kcp, &he, rU);
            }
        }
                
        complex_double_array_copy (tU, mU, NSB*NSB);
        
        /* find phase part using M */
        for (m=0;m<NSB;m++) {
            phi[m]=gsl_spline_eval (sM[m], kc, aTS)/epsilon;
        }
        
        for (m=0;m<NSB;m++) {
            complex double phase = cexp(-I*(phi[m]));
            for (n=0;n<NSB;n++) {
                tU[m*NSB+n]*=phase;
            }
        }
        
        complex_double_array_copy (&Lsb[NSB*NSB*q],tU,NSB*NSB);
        
        q--;
    }
        
    m_free(mU);
    m_free(tU);
    
    d_free(TSraw);
    
    for (m=0;m<2*NSB*NSB;m++) {
        gsl_interp_free (TS[m]);
    }
        
    for (m=0;m<NSB;m++) {
        gsl_spline_free(sM[m]);
    }
    
    gsl_odeiv_evolve_free (ee);
    gsl_odeiv_control_free (ce);
    gsl_odeiv_step_free (se);
    
    free(k);
    
    m_free(U);
    
    g_timer_stop(ztimer);
    
    return Lsb;
}

/* read out total time spent doing all this */
double zener_get_elapsed_time() {
    if (ztimer) {
        return g_timer_elapsed(ztimer, NULL);
    }
    else return 0.0;
}
