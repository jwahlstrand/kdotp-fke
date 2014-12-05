#include <complex.h>
#include <memory.h>
#include <glib.h>
#include <fftw3.h>
#include "three-vector.h"

#ifndef KDOTP_H
#define KDOTP_H

#define R 3.80998 // hbar^2/2m in eV A^2

#define CAL 0
#define OFFDIAG 1
#define INTRA 0
#define UCBANDS 0

#define m_free(x) fftw_free(x)
#define d_free(x) fftw_free(x)

static inline
void double_array_copy (double *array, const double * a, size_t n)
{
    memcpy (array, a, n*sizeof(double));
}

static inline double * double_array_calloc (size_t n) {
    double * array = (double *) fftw_malloc(sizeof(double) * n);
    memset(array, 0, sizeof(double)*n);
    return array;
}

static inline
double * double_array_clone (const double * a, size_t n)
{
    double * array = (double *) fftw_malloc(sizeof(double) * n);
    memcpy (array, a, n*sizeof(double));
    return array;
}

complex double  * complex_double_array_calloc (size_t n);

static inline
void complex_double_array_copy (complex double *array, const complex double * a, size_t n)
{
    memcpy (array, a, n*sizeof(complex double));
}

static inline
complex double * complex_double_array_clone (const complex double * a, size_t n)
{
    complex double * array = (complex double *) fftw_malloc(sizeof(complex double) * n);
    memcpy (array, a, n*sizeof(complex double));
    return array;
}

static inline
void complex_array_transpose (complex double * a, size_t N) {
    complex double * b = complex_double_array_clone (a, N*N);
    size_t i,j;
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++) {
            a[i+j*N]=b[j+i*N];
        }
    }
    m_free(b);
}

static inline void complex_array_zero (complex double * array , size_t n) {
    memset(array, 0, sizeof(complex double)*n);
}

void add_complex_double_arrays (complex double * darray, complex double * sarray, size_t n);

void scale_then_add_complex_double_arrays (complex double * darray, complex double * sarray, double * factor, size_t n);

void matrix_transform (complex double * B, const complex double * U, const complex double * A, size_t N);

void matrix_transform3 (complex double * Bx, complex double * By, complex double * Bz, const complex double * U, const complex double * Ax, const complex double * Ay, const complex double * Az, size_t N);

typedef struct {
	double P0, Q, P0p;
	double Eg, E0p;
	double D0, D0p;
	double G1L, G2L, G3L;
	double F;
	double db, Ck;
	double D1, D2, D3;
	double E6v,E6u,E3d,E5d,E6q;
	double P0_30,Qc,P3,P2,Ps,Pd,Qcd,P3d,P2d,Pu,Pm1;
} Material;

typedef struct {
    size_t total_bands, valence_bands, conduction_bands, first_valence, first_conduction;
    complex double * Px;
    complex double * Py;
    complex double * Pz;
    
    complex double * buffer1;
    
    complex double * (*H) (const ThreeVector * k);
    void (*dHdx) (complex double *, const ThreeVector * k);
    void (*dHdy) (complex double *, const ThreeVector * k);
    void (*dHdz) (complex double *, const ThreeVector * k);
    
    complex double * d2Hdx2;
    complex double * d2Hdy2;
    complex double * d2Hdz2;
    complex double * d2Hdxdy;
    complex double * d2Hdxdz;
    complex double * d2Hdydz;
    
    GList * bundles;
} Model;

typedef struct {
    GArray * list;
    size_t first;
} Bundle;

void k_dot_p_init ();
void k_dot_p_cleanup ();

GOptionGroup * k_dot_p_get_option_group ();

Material * material_new_by_name (char * name);

Model * model_new_by_name (char * name);

int model_has_remote_params ();
void model_free (Model * model);

#endif
