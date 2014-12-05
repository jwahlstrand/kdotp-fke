#include "zincblende.h"

extern Material * material;

static
complex double * Px2 ()
{   
	double P0 = material->P0;

    complex double * px = complex_double_array_calloc (2*2);
    
    px[0+1*2]=P0;
    px[1+0*2]=P0;

    return px;
}

static
complex double * Py2 ()
{    
	double P0 = material->P0;

    complex double * py = complex_double_array_calloc (2*2);
    
    py[0+1*2]=P0*I;
    py[1+0*2]=-P0*I;
    
    return py;
}

static
complex double * Pz2 ()
{   
	double P0 = material->P0;

    complex double * pz = complex_double_array_calloc (2*2);
    
    pz[0+1*2]=P0;
    pz[1+0*2]=P0;
    
    return pz;
}

complex double * H2pba (const ThreeVector * k)
{
    complex double * h = complex_double_array_calloc (2*2);
    
    double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
	double k2 = kx*kx+ky*ky+kz*kz;
    
    /*complex<double> * px = Px2();
    scale_then_add_complex_double_arrays (h, px, &kx, 4);
    fftw_free(px);
    complex<double> * py = Py2();
    scale_then_add_complex_double_arrays (h, py, &ky, 4);
    fftw_free(py);
    complex<double> * pz = Pz2();
    scale_then_add_complex_double_arrays (h, pz, &kz, 4);
    fftw_free(pz);*/
    
    h[0]=-(material->Eg)-R/0.45*k2;
    h[3]=R/0.08*k2;
    
    return h;
}

complex double * dHdx2pba (const ThreeVector * k)
{
    complex double * dh = complex_double_array_calloc (2*2);
    
    double kx = k->x[0];
    
    complex double * px = Px2();
    add_complex_double_arrays (dh, px, 2*2);
    fftw_free(px);
    
    dh[0]=-2*R*kx/0.45;
    dh[3]=2*R*kx/0.08;
    
    return dh;
}

complex double * dHdy2pba (const ThreeVector * k)
{
    complex double * dh = complex_double_array_calloc (2*2);
    
	double ky = k->x[1];
    
    complex double * py = Py2();
    add_complex_double_arrays (dh, py, 2*2);
    fftw_free(py);
    
    dh[0]=-2*R*ky/0.45;
    dh[3]=2*R*ky/0.08;
    
    return dh;
}

complex double * dHdz2pba (const ThreeVector * k)
{
    complex double * dh = complex_double_array_calloc (2*2);
    
	double kz = k->x[2];
    
    complex double * pz = Pz2();
    add_complex_double_arrays (dh, pz, 2*2);
    fftw_free(pz);
    
    dh[0]=-2*R*kz/0.45;
    dh[3]=2*R*kz/0.08;
    
    return dh;
}

Model * zincblende2pba_new ()
{
	Model * zincblende = g_slice_new(Model);
	zincblende->total_bands = 2;
	zincblende->valence_bands = 1;
	zincblende->conduction_bands = 1;
	zincblende->first_valence = 0;
	zincblende->first_conduction = 1;
	
	zincblende->Px = Px2 ();
	zincblende->Py = Py2 ();
	zincblende->Pz = Pz2 ();
	
	zincblende->H = H2pba;
	zincblende->dHdx = dHdx2pba;
	zincblende->dHdy = dHdy2pba;
	zincblende->dHdz = dHdz2pba;
	
	int n;
	GList * b = NULL;
	GArray * array = g_array_sized_new (FALSE, TRUE, sizeof(int), 1);
	n=0;
	g_array_append_val (array, n);  // not numbered in order of increasing energy!
	b=g_list_append (b, array);
	
	GArray * array4 = g_array_sized_new (FALSE, TRUE, sizeof(int), 1);
	n=1;
	g_array_append_val (array4, n);
	b=g_list_append (b, array4);
	
	zincblende->bundles = b;
	
	return zincblende;
}
