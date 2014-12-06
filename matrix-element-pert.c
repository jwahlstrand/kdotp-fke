#include <math.h>
#include <glib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include "matrix-element-pert.h"
#include "diagonalize.h"
#include "gsl_util.h"

/* calculate the matrix elements using perturbation theory */
/* we basically do step 1 here */

#define COARSENESS 4

#define DEBUG 0

extern Model * model;

typedef struct {
	ThreeVector * kperp;
	ThreeVector * dir;
} odeparams;

#define Fx(m,n) fc[(m)+(n)*NB]
#define Fy(m,n) fc[(m)+(n)*NB+N2]
#define Fz(m,n) fc[(m)+(n)*NB+2*N2]

#define O(m)    w[12*N2+m]
#define Op(m)   f[12*N2+m]

#define Wx(m,n) wc[(m)+(n)*NB]
#define Wy(m,n) wc[(m)+(n)*NB+N2]
#define Wz(m,n) wc[(m)+(n)*NB+2*N2]

#define Rxz(m,n) wc[(m)+(n)*NB+3*N2]
#define Ryz(m,n) wc[(m)+(n)*NB+4*N2]
#define Rzz(m,n) wc[(m)+(n)*NB+5*N2]

#define Rxzp(m,n) fc[(m)+(n)*NB+3*N2]
#define Ryzp(m,n) fc[(m)+(n)*NB+4*N2]
#define Rzzp(m,n) fc[(m)+(n)*NB+5*N2]

#if DEBUG
static
void degeneracy_alert (double kc, double f[], const double w[], int m, int n)
{
    complex double * fc = (complex double*) f;
    complex double * wc = (complex double*) w;
    int NB=model->total_bands;
    int N2=NB*NB;
    
    complex double * mat = complex_double_array_calloc (4);
    mat[0]=Wz(m,m);
    mat[1]=Wz(m,n);
    mat[2]=Wz(n,m);
    mat[3]=Wz(n,n);
    
    double * split = eigenvalues_herm (mat, 2);
    
    if (abs(Fz(m,n))>1e-9) {
        printf("%e %e -> %e %e\n",split[0],split[1], O(m), O(n));
        printf("%e %d %d: \t %e+%eI %e+%eI \n",kc, m, n, creal(Wz(m,m)), cimag(Wz(m,m)), creal(Wz(m,n)),cimag(Wz(m,n)));
        printf("%e \t \t %e+%eI %e+%eI \n", Op(m), creal(Wz(n,m)), cimag(Wz(n,m)), creal(Wz(n,n)),cimag(Wz(n,n)));
        printf("\n");
    }
}
#endif

/* matrix elements */

MatrixElement * matrix_element_new (const ThreeVector *k)
{
    int N = model->total_bands;
    MatrixElement * me = g_slice_new (MatrixElement);
    me->k = *k;
    me->energies = double_array_calloc (N);
    me->Wx = complex_double_array_calloc (N*N);
    me->Wy = complex_double_array_calloc (N*N);
    me->Wz = complex_double_array_calloc (N*N);
    return me;
}

void matrix_element_free (MatrixElement *me)
{
    d_free(me->energies);
    m_free(me->Wx);
    m_free(me->Wy);
    m_free(me->Wz);
    g_slice_free(MatrixElement, me);
    me = NULL;
}

void matrix_element_list_free (GList *list)
{
    GList *iter = list;
    while (iter) {
        matrix_element_free ( (MatrixElement *) iter->data);
        iter = g_list_next (iter);
    }
    g_list_free (list);
}

/* we could save memory in the list if we only save the info from the bands we will sum over */

MatrixElement * matrix_element_create (const ThreeVector *k, complex double * wallc)
{
    MatrixElement *me = matrix_element_new (k);
    me->kc=0;
    size_t NB = (model->total_bands);
    size_t N2 = (model->total_bands)*(model->total_bands);
    double * wall = (double *) wallc;
    
    complex_double_array_copy (me->Wx, wallc, N2);
    complex_double_array_copy (me->Wy, wallc+N2, N2);
    complex_double_array_copy (me->Wz, wallc+2*N2, N2);
    
    double_array_copy (me->energies, wall+12*N2, NB);
    return me;
}
