#include <complex.h>
#include <gsl/gsl_matrix.h>
#include "k-dot-p.h"
#include "three-vector.h"

#define KCMAX 0.5  // the highest value of KC that we calculate
#define KCMAXP 1.0 // the highest value of KC that we pad out to (for better resolution in the spectrum)

typedef struct {
    ThreeVector k;
    double kc;
    double *energies;
    complex double *Wx;
    complex double *Wy;
    complex double *Wz;
} MatrixElement;

MatrixElement * matrix_element_new (const ThreeVector *k);
void matrix_element_free (MatrixElement *me);
MatrixElement * matrix_element_create (const ThreeVector *k, complex double * wallc);
void matrix_element_list_free (GList *list);

MatrixElement * matrix_element_create_from_coeffs (const ThreeVector *k, complex double * callc, double kc);

/* in coeffs_pert.c */
GList * calc_w_phi_coeffs (const ThreeVector * kperp, const ThreeVector * normal);
void set_initial_condition2 (double * wall, gsl_matrix_complex * pxproj, gsl_matrix_complex * pyproj, gsl_matrix_complex * pzproj, complex double * rxzproj, complex double * ryzproj, complex double * rzzproj);

/* in coeffs_pert_other.c */
GList * calc_w_phi_constP (ThreeVector * kperp);
GList * calc_w_phi_KKP (ThreeVector * kperp, ThreeVector * dir);

double coeffs_pert_get_elapsed_time();
