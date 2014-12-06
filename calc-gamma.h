#include <glib.h>
#include "k-dot-p.h"
#include "three-vector.h"
#include <complex.h>

#ifndef _INC_CALC_GAMMA_H
#define _INC_CALC_GAMMA_H

/* represents Gamma for a particular pair of bands, and one kperp */
typedef struct {
        complex double * x;
        complex double * y;
        complex double * z;
        size_t Nkc;
        size_t max_N;  // maximum N that needs to be kept
        size_t min_N;
} Gamma;

Gamma * gamma_new(int Nkc);
void gamma_free(Gamma *g);

static inline
int near_degenerate (size_t m, size_t n)
{    
    if (labs(m-n)>1)
        return FALSE;
    if (m==n)
        return TRUE;
    if ((m % 2)==0 && n>m) {
        return TRUE;
    }
    if ((n % 2)==0 && m>n) {
        return TRUE;
    }
    return FALSE;
}

Bundle * find_bundle(size_t m);

size_t Nkc_from_epsilon (const double epsilon);

/* once we take the Fourier transform, we call it GammaHat */
typedef Gamma GammaHat;

/* return a list of GammaHats, one for each pair of bands, at a particular value of epsilon */
GammaHat ** calculate_gammahats (GList * matrix_elements, double epsilon, ThreeVector * dir, size_t max_N);
double calc_gamma_get_elapsed_time();

#endif
