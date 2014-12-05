#include "zincblende.h"
#include <gsl/gsl_math.h>

/* spherical 8 band model 8 */

#define epsilon 0.0 // was 1e-10 before
#define setmat(mat,i,j,val) mat[(j)+(i)*8]=(val)
#define getmat(mat,i,j) mat[(j)+(i)*8]

extern Material * material;
extern Model * model;

static
complex double * H8sph (const ThreeVector * k)
{
	complex double * h = complex_double_array_calloc (8*8);
	
	double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
	double k2 = kx*kx+ky*ky+kz*kz;  // calculate these now
    double Ek=R*k2;
    double Ez=R*kz*kz;
    double E2zmxy=R*((kz+kx)*(kz-kx)+(kz+ky)*(kz-ky));
    
    double EP=(material->P0)*(material->P0)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg));
    double g2=(material->G2L)-EP/(6*(material->Eg));
    double g3=(material->G3L)-EP/(6*(material->Eg));
    
    double g23=(2*g2+3*g3)/5;
    
    scale_then_add_complex_double_arrays (h, model->Px, &kx, 64);
    scale_then_add_complex_double_arrays (h, model->Py, &ky, 64);
    scale_then_add_complex_double_arrays (h, model->Pz, &kz, 64);
	
	/* diagonal elements */
    
    complex double a = Ek*2*(material->F)+Ek;  // conduction band
    setmat(h, 0, 0, a);
    setmat(h, 4, 4, a+epsilon);
    
    a=-(material->Eg)-Ek*(g1-g23)-3*Ez*g23;
    setmat(h, 1, 1, a);
    setmat(h, 5, 5, a+epsilon);
    
    a=-(material->Eg)-Ek*(g1+g23)+3*Ez*g23;
    setmat(h, 2, 2, a);
    setmat(h, 6, 6, a+epsilon);
    
    a=-((material->Eg)+(material->D0))-Ek*g1;
    setmat(h, 3, 3, a);
    setmat(h, 7, 7, a+epsilon);
    
    complex double kp = (kx + I*ky)/M_SQRT2;
    complex double km=conj(kp);
    
    setmat(h, 1, 2, M_SQRT3*R*(kx+ky)*(kx-ky)*g23-I*M_SQRT3*R*kx*ky*2.0*g23);
    setmat(h, 1, 3, M_SQRT2*g23*E2zmxy);
    setmat(h, 1, 6, 2*M_SQRT3*M_SQRT2*g23*R*kz*kp);
    setmat(h, 1, 7, 6*g23*R*kz*km);
    
    setmat(h, 2, 1, conj(getmat(h,1,2)));
    setmat(h, 2, 3, M_SQRT3*M_SQRT2*R*(kx+ky)*(kx-ky)*g23+I*M_SQRT3*M_SQRT2*R*kx*ky*2.0*g23);
    setmat(h, 2, 5, 2*M_SQRT3*M_SQRT2*g23*R*kz*kp);
    setmat(h, 2, 7, -2*M_SQRT3*g23*R*kz*kp);
    
    setmat (h, 3, 1, conj(getmat(h, 1, 3)));
    setmat (h, 3, 2, conj(getmat(h, 2, 3)));
    setmat (h, 3, 5, -6.0*g23*R*kz*km);
    setmat (h, 3, 6, getmat(h,2,7));
    
    setmat (h, 5, 2, conj(getmat(h, 2, 5)));
    setmat (h, 5, 3, conj(getmat(h, 3, 5)));
    setmat (h, 5, 6, M_SQRT3*R*(kx+ky)*(ky-kx)*g23-I*M_SQRT3*R*kx*ky*2.0*g23);
    setmat (h, 5, 7, getmat(h,1,3));
    
    setmat (h, 6, 1, conj(getmat(h, 1, 6)));
    setmat (h, 6, 3, conj(getmat(h, 3, 6)));
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    setmat (h, 6, 7, M_SQRT3*M_SQRT2*R*(kx+ky)*(ky-kx)*g23+I*M_SQRT3*M_SQRT2*R*kx*ky*2.0*g23);
    
    setmat (h, 7, 1, conj(getmat(h, 1, 7)));
    setmat (h, 7, 2, conj(getmat(h, 2, 7)));
    setmat (h, 7, 5, conj(getmat(h, 5, 7)));
    setmat (h, 7, 6, conj(getmat(h, 6, 7)));

	return h;
}

static
complex double * d2Hdx2_8sph (){
    complex double * h = complex_double_array_calloc (8*8);
    
    double EP=(material->P0)*(material->P0)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg));
    double g2=(material->G2L)-EP/(6*(material->Eg));
    double g3=(material->G3L)-EP/(6*(material->Eg));
    
    double g23=(2*g2+3*g3)/5;
    
    complex double a = 2*R*2*(material->F)+2*R;  // conduction band
    setmat(h, 0, 0, a);
    setmat(h, 4, 4, a);
    
    a=-2*R*(g1-g23);
    setmat(h, 1, 1, a);
    setmat(h, 5, 5, a);
    
    a=-2*R*(g1+g23);
    setmat(h, 2, 2, a);
    setmat(h, 6, 6, a);
    
    a=-2*R*g1;
    setmat(h, 3, 3, a);
    setmat(h, 7, 7, a);
        
    setmat(h, 1, 2, M_SQRT3*2*R*g23);
    setmat(h, 1, 3, -2*R*M_SQRT2*g23);
    
    setmat(h, 2, 1, conj(getmat(h,1,2)));
    setmat(h, 2, 3, M_SQRT3*M_SQRT2*2*R*g23);
    
    setmat (h, 3, 1, conj(getmat(h, 1, 3)));
    setmat (h, 3, 2, conj(getmat(h, 2, 3)));
    
    setmat (h, 5, 6, -2*M_SQRT3*R*g23);
    setmat (h, 5, 7, getmat(h,1,3));
    
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    setmat (h, 6, 7, -2*M_SQRT3*M_SQRT2*R*g23);
    
    setmat (h, 7, 5, conj(getmat(h, 5, 7)));
    setmat (h, 7, 6, conj(getmat(h, 6, 7)));
    
    return h;
}

static
complex double * d2Hdy2_8sph (){
    complex double * h = complex_double_array_calloc (8*8);
    
    double EP=(material->P0)*(material->P0)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg));
    double g2=(material->G2L)-EP/(6*(material->Eg));
    double g3=(material->G3L)-EP/(6*(material->Eg));
    
    double g23=(2*g2+3*g3)/5;
    
    complex double a = 2*R*2*(material->F)+2*R;  // conduction band
    setmat(h, 0, 0, a);
    setmat(h, 4, 4, a);
    
    a=-2*R*(g1-g23);
    setmat(h, 1, 1, a);
    setmat(h, 5, 5, a);
    
    a=-2*R*(g1+g23);
    setmat(h, 2, 2, a);
    setmat(h, 6, 6, a);
    
    a=-2*R*g1;
    setmat(h, 3, 3, a);
    setmat(h, 7, 7, a);
    
    setmat(h, 1, 2, -2*M_SQRT3*R*g23);
    setmat(h, 1, 3, -2*M_SQRT2*g23*R);
    
    setmat(h, 2, 1, conj(getmat(h,1,2)));
    setmat(h, 2, 3, -2*M_SQRT3*M_SQRT2*R*g23);
    
    setmat (h, 3, 1, conj(getmat(h, 1, 3)));
    setmat (h, 3, 2, conj(getmat(h, 2, 3)));
    
    setmat (h, 5, 6, 2*M_SQRT3*R*g23);
    setmat (h, 5, 7, getmat(h,1,3));
    
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    setmat (h, 6, 7, 2*M_SQRT3*M_SQRT2*R*g23);
    
    setmat (h, 7, 5, conj(getmat(h, 5, 7)));
    setmat (h, 7, 6, conj(getmat(h, 6, 7)));
    
    return h;
}

static
complex double * d2Hdz2_8sph (){
    complex double * h = complex_double_array_calloc (8*8);
    
    double EP=(material->P0)*(material->P0)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg));
    double g2=(material->G2L)-EP/(6*(material->Eg));
    double g3=(material->G3L)-EP/(6*(material->Eg));
    
    double g23=(2*g2+3*g3)/5;
    
    complex double a = 2*R*2*(material->F)+2*R;  // conduction band
    setmat(h, 0, 0, a);
    setmat(h, 4, 4, a);
    
    a=-2*R*(g1-g23)-6*R*g23;
    setmat(h, 1, 1, a);
    setmat(h, 5, 5, a);
    
    a=-2*R*(g1+g23)+6*R*g23;
    setmat(h, 2, 2, a);
    setmat(h, 6, 6, a);
    
    a=-2*R*g1;
    setmat(h, 3, 3, a);
    setmat(h, 7, 7, a);
    
    setmat(h, 1, 3, 4*R*M_SQRT2*g23);
        
    setmat (h, 3, 1, conj(getmat(h, 1, 3)));
    
    setmat (h, 5, 7, getmat(h,1,3));
        
    setmat (h, 7, 5, conj(getmat(h, 5, 7)));
    
    return h;
}

static
complex double * d2Hdxdy_8sph (){
    complex double * h = complex_double_array_calloc (8*8);
    
    double EP=(material->P0)*(material->P0)/R;
    double g2=(material->G2L)-EP/(6*(material->Eg));
    double g3=(material->G3L)-EP/(6*(material->Eg));
    
    double g23=(2*g2+3*g3)/5;
    
    setmat(h, 1, 2, -I*M_SQRT3*R*2.0*g23);
    
    setmat(h, 2, 1, conj(getmat(h,1,2)));
    setmat(h, 2, 3, I*M_SQRT3*M_SQRT2*R*2.0*g23);
    
    setmat (h, 3, 2, conj(getmat(h, 2, 3)));
    
    setmat (h, 5, 6, -I*M_SQRT3*R*2.0*g23);
    
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    setmat (h, 6, 7, I*M_SQRT3*M_SQRT2*R*2.0*g23);
    
    setmat (h, 7, 6, conj(getmat(h, 6, 7)));
    
    return h;
}

static
complex double * d2Hdxdz_8sph (){
    complex double * h = complex_double_array_calloc (8*8);
    
    double EP=(material->P0)*(material->P0)/R;
    double g2=(material->G2L)-EP/(6*(material->Eg));
    double g3=(material->G3L)-EP/(6*(material->Eg));
    
    double g23=(2*g2+3*g3)/5;
    
    setmat(h, 1, 6, 2*M_SQRT3*g23*R);
    setmat(h, 1, 7, 6*g23*R/M_SQRT2);
    
    setmat(h, 2, 5, 2*M_SQRT3*g23*R);
    setmat(h, 2, 7, -2*M_SQRT3/M_SQRT2*g23*R);
    
    setmat (h, 3, 5, -6.0*g23*R/M_SQRT2);
    setmat (h, 3, 6, getmat(h,2,7));
    
    setmat (h, 5, 2, conj(getmat(h, 2, 5)));
    setmat (h, 5, 3, conj(getmat(h, 3, 5)));
    
    setmat (h, 6, 1, conj(getmat(h, 1, 6)));
    setmat (h, 6, 3, conj(getmat(h, 3, 6)));
    
    setmat (h, 7, 1, conj(getmat(h, 1, 7)));
    setmat (h, 7, 2, conj(getmat(h, 2, 7)));
    
    return h;
}

static
complex double * d2Hdydz_8sph (){
    complex double * h = complex_double_array_calloc (8*8);
    
    double EP=(material->P0)*(material->P0)/R;
    double g2=(material->G2L)-EP/(6*(material->Eg));
    double g3=(material->G3L)-EP/(6*(material->Eg));
    
    double g23=(2*g2+3*g3)/5;
    
    setmat(h, 1, 6, 2*M_SQRT3*g23*R*I);
    setmat(h, 1, 7, 6*g23*R*(-I)/M_SQRT2);
    
    setmat(h, 2, 5, 2*M_SQRT3*g23*R*I);
    setmat(h, 2, 7, -2*M_SQRT3*g23*R*I/M_SQRT2);
    
    setmat (h, 3, 5, -6.0*g23*R*(-I)/M_SQRT2);
    setmat (h, 3, 6, getmat(h,2,7));
    
    setmat (h, 5, 2, conj(getmat(h, 2, 5)));
    setmat (h, 5, 3, conj(getmat(h, 3, 5)));
    
    setmat (h, 6, 1, conj(getmat(h, 1, 6)));
    setmat (h, 6, 3, conj(getmat(h, 3, 6)));
    
    setmat (h, 7, 1, conj(getmat(h, 1, 7)));
    setmat (h, 7, 2, conj(getmat(h, 2, 7)));
    
    return h;
}

static
complex double * dHdx8 (const ThreeVector * k)
{
    complex double * h = complex_double_array_calloc (8*8);
    
    add_complex_double_arrays (h, model->Px, 64);
    
    double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
    scale_then_add_complex_double_arrays (h, model->d2Hdx2, &kx, 64);
    scale_then_add_complex_double_arrays (h, model->d2Hdxdy, &ky, 64);
    scale_then_add_complex_double_arrays (h, model->d2Hdxdz, &kz, 64);
    
    return h;
}

static
complex double * dHdy8 (const ThreeVector * k)
{
    complex double * h = complex_double_array_calloc (8*8);
    
    add_complex_double_arrays (h, model->Py, 64);
    
    double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
    scale_then_add_complex_double_arrays (h, model->d2Hdxdy, &kx, 64);
    scale_then_add_complex_double_arrays (h, model->d2Hdy2, &ky, 64);
    scale_then_add_complex_double_arrays (h, model->d2Hdydz, &kz, 64);
    
    return h;
}

static
complex double * dHdz8 (const ThreeVector * k)
{
    complex double * h = complex_double_array_calloc (8*8);
    
    add_complex_double_arrays (h, model->Pz, 64);
    
    double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
    scale_then_add_complex_double_arrays (h, model->d2Hdxdz, &kx, 64);
    scale_then_add_complex_double_arrays (h, model->d2Hdydz, &ky, 64);
    scale_then_add_complex_double_arrays (h, model->d2Hdz2, &kz, 64);
    
    return h;
}

Model * zincblende8sph_new ()
{
	Model * zincblende = g_slice_new(Model);
	zincblende->total_bands = 8;
	zincblende->valence_bands = 6;
	zincblende->conduction_bands = 2;
	zincblende->first_valence = 0;
	zincblende->first_conduction = 6;
	
	zincblende->Px = Px8 ();
	zincblende->Py = Py8 ();
	zincblende->Pz = Pz8 ();
	
	zincblende->H = H8sph;
	
	zincblende->dHdx = dHdx8;
	zincblende->dHdy = dHdy8;
	zincblende->dHdz = dHdz8;
	
	zincblende->d2Hdx2 = d2Hdx2_8sph ();
	zincblende->d2Hdy2 = d2Hdy2_8sph ();
	zincblende->d2Hdz2 = d2Hdz2_8sph ();
	zincblende->d2Hdxdy = d2Hdxdy_8sph ();
	zincblende->d2Hdxdz = d2Hdxdz_8sph ();
	zincblende->d2Hdydz = d2Hdydz_8sph ();
	    
    bundle8(zincblende);
	
	return zincblende;
}
