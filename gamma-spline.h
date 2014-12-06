#include "config.h"
#include <glib.h>
#include <gsl/gsl_spline.h>
#include "calc-gamma.h"

#ifndef GAMMA_SPLINE_H
#define GAMMA_SPLINE_H

typedef struct {
        double *k, *phi, *Wxr, *Wxi, *Wyr, *Wyi, *Wzr, *Wzi;
        
        gsl_interp * phi_spline;
        gsl_interp_accel * phi_accel;
        gsl_interp * wx_spline_re;
        gsl_interp_accel * w_accel;
        gsl_interp * wx_spline_im;
    
        gsl_interp * wy_spline_re;
        gsl_interp * wy_spline_im;
    
        gsl_interp * wz_spline_re;
        gsl_interp * wz_spline_im;
} GammaSpline;

GammaSpline * gamma_spline_new(int size);
void gamma_spline_reset (GammaSpline *gs, GList * const sd, int mp, int np, gboolean sort);
double * gamma_spline_energy(GammaSpline *gs, GList * const sd, double kc);
complex double * gamma_spline_V(GammaSpline *gs, GList * const sd, double kc, int component);
Gamma * gamma_spline_calc_little_gamma (GammaSpline *gs, double * omega);
Gamma * gamma_spline_calc_gamma (const GammaSpline *gs, double epsilon);

void gamma_spline_free(GammaSpline *gs);

double gamma_spline_get_elapsed_time();

#endif
