#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include "matrix-element-pert.h"
#include "gamma-spline.h"
#include "calc-theta.h"
#include "absorption.h"

#define DEBUG 0

extern Model * model;
extern Material * material;

/* Spectrum:  base class of OnePhoton, QIC, and TwoPhoton */

void spectrum_set_N1_N2_from_epsilon (OnePhoton *op, double epsilon, int factor) {
	int q, length;
	double ostep;
	if (fabs(epsilon)<1e-6) {
	    length = 512;
	}
	else {
    	length = Nkc_from_epsilon(epsilon);
    }
    //g_message("length is %d",length);
    double min_energy = (3.0/4.0*material->Eg)/factor;
    double max_energy = (material->Eg*1.5)/factor;
    if (fabs(epsilon)<1e-6) {
        ostep=(max_energy-min_energy)/length;
    }
    else {
        double kstep=4*KCMAX/length;
        ostep=fabs(epsilon)/kstep*1240.7/3e5/length;
    }
    //g_message("ostep is %e",ostep);
	op->N1 = (int) floor(min_energy/ostep);
	op->N2 = (int) ceil(max_energy/ostep);
	const int N = op->N2 - op->N1;
	op->o=double_array_calloc(N);
    for (q=0;q<N;q++) {
        op->o[q]=((double)q+op->N1)*ostep;  // only include energies between min_energy and max_energy
    }
  	//if (factor==2)
    //	printf("N1 is %d, N2 is %d, o[0] is %e\n",N1,N2,o[0]);
    op->e=epsilon;
}

/* OnePhoton */

/* Represents a one-photon absorption spectrum, including off-diagonal components */

OnePhoton * one_photon_new (int N)
{
    int q;
    OnePhoton *op = g_slice_new(OnePhoton);
    double min_energy = 1.0;
    double max_energy = 2.0;
    op->e=0;
	double ostep=(max_energy-min_energy)/N;
	op->N1 = 0;
	op->N2 = N;
	op->o=double_array_calloc(N);
    for (q=0;q<N;q++) {
        op->o[q]=min_energy+((double)q)*ostep;  // only include energies between min_energy and max_energy
    }
    op->data = (Tensor2 *) malloc (sizeof(Tensor2)*N);
    memset(op->data, 0, sizeof(Tensor2)*N);
    return op;
}

OnePhoton * one_photon_epsilon (double epsilon)
{
    OnePhoton *op = g_slice_new(OnePhoton);

    spectrum_set_N1_N2_from_epsilon(op, epsilon, 1);
	const int N = op->N2 - op->N1;
    op->data = (Tensor2 *) malloc (sizeof(Tensor2)*N);
    memset(op->data, 0, sizeof(Tensor2)*N);
    return op;
}

void one_photon_free (OnePhoton *op)
{
    free(op->data);
    g_slice_free(OnePhoton,op);
}

/* scale each component by weight */

void one_photon_scale (OnePhoton *op, double weight)
{
    int q,p;
    int N = op->N2 - op->N1;
    for (q=0;q<N;q++) {
        for (p=0;p<3;p++) {
            op->data[q].diag[p]=op->data[q].diag[p]*weight;
            op->data[q].offdiag[p]=op->data[q].offdiag[p]*weight;
        }
    }
}

/* scale another spectrum by weight and add to this spectrum */

void one_photon_add_and_scale (OnePhoton *op1, OnePhoton *op2, double weight)
{
    int q,p;
    int N = op1->N2 - op1->N1;
    for (q=0;q<N;q++) {
        for (p=0;p<3;p++) {
            op1->data[q].diag[p]+=op2->data[q].diag[p]*weight;
            op1->data[q].offdiag[p]+=op2->data[q].offdiag[p]*weight;
        }
    }
}

/* write spectrum to a binary file */
/* we write only doubles to keep decoding simple */

void one_photon_write (OnePhoton *op, FILE* file)
{
    double t = 1.0; //signifies that this is one-photon absorption
    double N = (double) (op->N2 - op->N1);
    size_t size = sizeof(double);
    fwrite(&t, size, 1, file);  // write type code as a double
    fwrite(&N, size, 1, file);  // write N as a double
    double omegastep = op->o[1]-op->o[0];
    double realo1;
    if (fabs(op->e)<1e-6)
        realo1 = op->o[0]+omegastep/2.0;  // we want energy to be the center of the bin for the zero field case
    else
        realo1 = op->o[0];
    fwrite(&realo1, size, 1, file);        // write omega1 as double
    fwrite(&omegastep, size, 1, file);   // write omegastep as double
    fwrite(&op->e, size, 1, file);           // write epsilon as double
    fwrite(&op->data[0], sizeof(Tensor2), (size_t) N, file);
}

/* This takes Gammahat and calculates a spectrum from it */

void one_photon_calc_absorption_component (OnePhoton *op, GammaHat * gamma) {
	int q;
	double omegastep = op->o[1]-op->o[0];
	for (q=op->N1;q<op->N2;q++) {
    	double omega = q*omegastep;
        op->data[q-op->N1].diag[0]+=gsl_pow_2(cabs(gamma->x[q]))/pow(omega,2);
        op->data[q-op->N1].diag[1]+=gsl_pow_2(cabs(gamma->y[q]))/pow(omega,2);
        op->data[q-op->N1].diag[2]+=gsl_pow_2(cabs(gamma->z[q]))/pow(omega,2);
        op->data[q-op->N1].offdiag[0]+=conj(gamma->x[q])*gamma->y[q]/pow(omega,2);
        op->data[q-op->N1].offdiag[1]+=conj(gamma->x[q])*gamma->z[q]/pow(omega,2);
        op->data[q-op->N1].offdiag[2]+=conj(gamma->y[q])*gamma->z[q]/pow(omega,2);
	}
}

void one_photon_incr_nofield_absorption_component (OnePhoton *op, GammaHat *gamma, double omega, int k) {
    double omegastep = op->o[1]-op->o[0];
    int q = floor(omega/omegastep);
    if (q>=op->N1 && q<op->N2) {
        op->data[q-op->N1].diag[0]+=gsl_pow_2(cabs(gamma->x[k]))/pow(omega,2);
        op->data[q-op->N1].diag[1]+=gsl_pow_2(cabs(gamma->y[k]))/pow(omega,2);
        op->data[q-op->N1].diag[2]+=gsl_pow_2(cabs(gamma->z[k]))/pow(omega,2);
        op->data[q-op->N1].offdiag[0]+=conj(gamma->x[k])*gamma->y[k]/pow(omega,2);
        op->data[q-op->N1].offdiag[1]+=conj(gamma->x[k])*gamma->z[k]/pow(omega,2);
        op->data[q-op->N1].offdiag[2]+=conj(gamma->y[k])*gamma->z[k]/pow(omega,2);
        
    }
}

/****************************/

static gboolean include_zero=TRUE;

static
double get_eps (double max_eps, int total, int p) {
	if (include_zero) {
    	if (total <= 1)
    		return 0.0;
    	double deps = max_eps/(total-1);
    	return p*deps;
    }
    else {
        if (total <= 1)
    		return max_eps;
    	double deps = max_eps/(total);
    	return (p+1)*deps;
    }
}


/* creates a spectrum of all zeros (for bad boxes) */

OnePhoton ** blank_abs (const double max_eps, const int num, gboolean incl_zero) {
    OnePhoton ** as = (OnePhoton **) malloc(sizeof(OnePhoton *)*num);
    
    include_zero = incl_zero;
  	
  	int p;
  	double epsilon;
  	for (p=0;p<num;p++) {
		epsilon = get_eps(max_eps, num, p);
		as[p]=one_photon_epsilon (epsilon);
	}
    return as;
}

/* take a list of matrix elements and calculate one-photon absorption spectra for a bunch of epsilons */

OnePhoton ** calculate_absorption (GList * list, const double max_eps, const int num, ThreeVector * dir, gboolean incl_zero) {
    int p, q, num_e;
    int Nkc;
    double en, starte, ende, de, epsilon;
    
    double * omega;
    
    include_zero = incl_zero;
    
  	OnePhoton ** as = (OnePhoton **) malloc(sizeof(OnePhoton *)*num);
    
	starte = 3.0/4.0*material->Eg;
	ende = 1.5*(material->Eg);
	
  	for (p=0;p<num;p++) {
		epsilon = get_eps(max_eps, num, p);
		as[p]=one_photon_epsilon (epsilon);
	}
	
	num_e = (as[0]->N2-as[0]->N1);
	de=(ende-starte)/num_e;
	
	int NB = model->total_bands;
	
	GammaHat ** gh = NULL;
    
    for (p=0;p<num;p++) {
    	epsilon = get_eps(max_eps, num, p);
    	
    	g_message("calculating for eps=%f",epsilon);
    	
    	if (fabs(epsilon)<1e-5) {
            size_t nv, nc;
            GammaSpline * gs = gamma_spline_new (g_list_length(list));
            
            for (nv=model->first_valence;nv<model->first_conduction;nv++) {
            //for (nv=0;nv<2;nv++) {
                for (nc=model->first_conduction;nc<model->total_bands;nc++) {
                    gamma_spline_reset(gs,list, nv, nc, TRUE);
                    
            	    omega = double_array_calloc (16384);
            	    Nkc=16384;
        		    Gamma * g = gamma_spline_calc_little_gamma (gs,omega); /* omega gets set here */
        		            		    
        		    for (q=0;q<Nkc;q++) {
                        en=omega[q]*1240.7/3e5/2/M_PI;
                        
                        one_photon_incr_nofield_absorption_component(as[0],g, en, q);
        				
                    }
                    gamma_free(g);
                    d_free(omega);
                }
            }
            gamma_spline_free(gs);
    	}
    	else {
    	    gh=calculate_gammahats (list, epsilon, dir, 10000);
			//printf("size is %d, Nkc is %d\n",as->one[p]->N2, Nkc_from_epsilon(epsilon));
            size_t v, c;
            for (v=(model->first_valence);v<(model->first_conduction);v++) {
            //for (v=0;v<2;v++) {
                GammaHat * g;
                for (c=(model->first_conduction);c<(model->first_conduction+model->conduction_bands);c++) {
                    g = gh[c+NB*v];
                    
    			    one_photon_calc_absorption_component(as[p],g);
    		        
    		    }
    		}
			
			for (q=0;q<NB*NB;q++) {
			    if (gh[q])
    			    gamma_free(gh[q]);
			}
			//delete[] gh;
        }
    }
    
    double dkc=2*KCMAX/16384;
    
    /* we calculate the first order susceptibility in SI units */
    /* we calculate the QUIC injection tensor in SI units */
    /* we calculate the second-order injection tensor in SI units */
    
    if (include_zero) {
        one_photon_scale(as[0],dkc/de*2.2918);  // last factor related to e^2, hbar, 2*pi, etc. (see units write-up)
    }
    int start = include_zero ? 1 : 0;
    for (p=start;p<num;p++) {
        epsilon = get_eps(max_eps, num, p);
        
    	one_photon_scale(as[p],fabs(epsilon)*0.00947802);
	}
        
    //gsl_vector_free(omega);
		
	return as;
}

