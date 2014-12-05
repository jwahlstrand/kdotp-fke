#include "zincblende.h"
#include <gsl/gsl_math.h>

#define epsilon 0.0//1e-10
#define setmat(mat,i,j,val) mat[(i)+(j)*14]=(val)
#define getmat(mat,i,j) mat[(i)+(j)*14]

extern Material * material;

static
complex double * H14nr (const ThreeVector * k)
{
	complex double * h = complex_double_array_calloc (14*14);
	
	double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
	double k2 = kx*kx+ky*ky+kz*kz;  // calculate these now
    double Ek=R*k2;
        
    double Eg = material->Eg;
    double E1 = material->E0p - material->Eg;
    
    double G0 = -Eg-material->D0;
    double G1 = E1+material->D0p;
    
    complex double * px = Px14();
    scale_then_add_complex_double_arrays (h, px, &kx, 196);
    fftw_free(px);
    
    complex double * py = Py14();
    scale_then_add_complex_double_arrays (h, py, &ky, 196);
    fftw_free(py);
    
    complex double * pz = Pz14();
    scale_then_add_complex_double_arrays (h, pz, &kz, 196);
    fftw_free(pz);
	
	/* diagonal elements */
	
	complex double a = G1+Ek;
    setmat (h, 0, 0, a);
    setmat (h, 1, 1, a+epsilon);
    setmat (h, 7, 7, a);
    setmat (h, 8, 8, a+epsilon);
    
    a = E1+Ek;
    setmat (h, 2, 2, a);
    setmat (h, 9, 9, a+epsilon);
            
    a = Ek;  // conduction band
    setmat(h, 3, 3, a);
    setmat(h, 10, 10, a+epsilon);
    
    a=-Eg+Ek;
    setmat(h, 4, 4, a);
    setmat(h, 11, 11, a+epsilon);
    
    a=-Eg+Ek;
    setmat(h, 5, 5, a);
    setmat(h, 12, 12, a+epsilon);
    
    a=G0+Ek;
    setmat(h, 6, 6, a);
    setmat(h, 13, 13, a+epsilon);
    
	return h;
}

static
complex double * d2Hdall2_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    complex double a = 2*R;
    setmat (h, 0, 0, a);
    setmat (h, 1, 1, a);
    setmat (h, 7, 7, a);
    setmat (h, 8, 8, a);
    setmat (h, 2, 2, a);
    setmat (h, 9, 9, a);
    
    a = 2*R;  // conduction band
    setmat(h, 3, 3, a);
    setmat(h, 10, 10, a);
    
    a=2*R;
    setmat(h, 4, 4, a);
    setmat(h, 11, 11, a);
    
    a=2*R;
    setmat(h, 5, 5, a);
    setmat(h, 12, 12, a);
    
    a=2*R;
    setmat(h, 6, 6, a);
    setmat(h, 13, 13, a);
    
    return h;
}

static
complex double * zero_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    return h;
}

static
complex double * dHdx14 (const ThreeVector * k)
{
    complex double * h = complex_double_array_calloc (196);
    
    double one = 1;
    
    complex double * px = Px14();
    scale_then_add_complex_double_arrays (h, px, &one, 196);
    fftw_free(px);
    
    double kx = k->x[0];
	
	px = d2Hdall2_14();
    scale_then_add_complex_double_arrays (h, px, &kx, 196);
    fftw_free(px);
    
    return h;
}

static
complex double * dHdy14 (const ThreeVector * k)
{
    complex double * h = complex_double_array_calloc (196);
    
    double one = 1;
    
    complex double * py = Py14();
    scale_then_add_complex_double_arrays (h, py, &one, 196);
    fftw_free(py);
    
	double ky = k->x[1];
	
    py = d2Hdall2_14();
    scale_then_add_complex_double_arrays (h, py, &ky, 196);
    fftw_free(py);
    
    return h;
}

static
complex double * dHdz14 (const ThreeVector * k)
{
    complex double * h = complex_double_array_calloc (196);
    
    double one = 1;
    
    complex double * pz = Pz14();
    scale_then_add_complex_double_arrays (h, pz, &one, 196);
    fftw_free(pz);
    
	double kz = k->x[2];
	
    pz = d2Hdall2_14();
    scale_then_add_complex_double_arrays (h, pz, &kz, 196);
    fftw_free(pz);
    
    return h;
}

Model * zincblende14nr_new ()
{
	Model * zincblende = g_slice_new(Model);
	zincblende->total_bands = 14;
	zincblende->valence_bands = 6;
	zincblende->conduction_bands = 2;
	zincblende->first_valence = 0;
	zincblende->first_conduction = 6;
	
	zincblende->Px = Px14 ();
	zincblende->Py = Py14 ();
	zincblende->Pz = Pz14 ();
	
	zincblende->H = H14nr;
	
	zincblende->dHdx = dHdx14;
	zincblende->dHdy = dHdy14;
	zincblende->dHdz = dHdz14;
	
	zincblende->d2Hdx2 = d2Hdall2_14 ();
	zincblende->d2Hdy2 = d2Hdall2_14 ();
	zincblende->d2Hdz2 = d2Hdall2_14 ();
	zincblende->d2Hdxdy = zero_14 ();
	zincblende->d2Hdxdz = zero_14 ();
	zincblende->d2Hdydz = zero_14 ();
	
    bundle14(zincblende);
	
	return zincblende;
}
