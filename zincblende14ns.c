#include "zincblende.h"
#include <gsl/gsl_math.h>

#define epsilon 0.0//1e-10
#define setmat(mat,i,j,val) mat[(j)+(i)*14]=(val)
#define getmat(mat,i,j) mat[(j)+(i)*14]

extern Material * material;

complex double * Px14 ()
{   
	double P0 = material->P0;
	double Q = material->Q;
	double P0p = material->P0p;

    complex double * px = complex_double_array_calloc (14*14);
    
    #include "px14.txt"

    return px;
}

complex double * Py14 ()
{    
	double P0 = material->P0;
	double Q = material->Q;
	double P0p = material->P0p;

    complex double * py = complex_double_array_calloc (14*14);
    
    #include "py14.txt"
    
    complex_array_transpose(py,14);
    
    return py;
}
    
complex double * Pz14 ()
{   
	double P0 = material->P0;
	double Q = material->Q;
	double P0p = material->P0p;

    complex double * pz = complex_double_array_calloc (14*14);
    
    #include "pz14.txt"
    
    return pz;
}

static
complex double * H14ns (const ThreeVector * k)
{
	complex double * h = complex_double_array_calloc (14*14);
	
	double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
	double k2 = kx*kx+ky*ky+kz*kz;  // calculate these now
    double Ek=R*k2;
    double Ez=R*kz*kz;
    double E2zmxy=R*((kz+kx)*(kz-kx)+(kz+ky)*(kz-ky));
    
    double EP=(material->P0)*(material->P0)/R;
    double EQ=(material->Q)*(material->Q)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg))-EQ/(3*(material->E0p))-EQ/(3*((material->E0p)+(material->D0p)));
    double g2=(material->G2L)-EP/(6*(material->Eg))+EQ/(6*(material->E0p));
    double g3=(material->G3L)-EP/(6*(material->Eg))-EQ/(6*(material->E0p));
    
    double Eg = material->Eg;
    double E1 = material->E0p - material->Eg;
    
    double G0 = -Eg-material->D0;
    double G1 = E1+material->D0p;
    
    /*printf("E0 is %f, ",-Eg);
    printf("E1 is %f, ",E1);
    printf("G0 is %f, ",G0);
    printf("G1 is %f\n", G1);*/
    
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
            
    a = Ek*2*(material->F)+Ek;  // conduction band
    setmat(h, 3, 3, a);
    setmat(h, 10, 10, a+epsilon);
    
    a=-Eg-Ek*(g1-g2)-3*Ez*g2;
    setmat(h, 4, 4, a);
    setmat(h, 11, 11, a+epsilon);
    
    a=-Eg-Ek*(g1+g2)+3*Ez*g2;
    setmat(h, 5, 5, a);
    setmat(h, 12, 12, a+epsilon);
    
    a=G0-Ek*g1;
    setmat(h, 6, 6, a);
    setmat(h, 13, 13, a+epsilon);
    
    complex double kp = (kx + I*ky)/M_SQRT2;
    complex double km=conj(kp);
    
    setmat(h, 4, 5, M_SQRT3*R*(kx+ky)*(kx-ky)*g2-I*M_SQRT3*R*kx*ky*2.0*g3);
    setmat(h, 4, 6, M_SQRT2*g2*E2zmxy);
    setmat(h, 4, 12, 2*M_SQRT3*M_SQRT2*g3*R*kz*kp);
    setmat(h, 4, 13, 6*g3*R*kz*km);
    
    setmat(h, 5, 4, conj(getmat(h,4,5)));
    setmat(h, 5, 6, M_SQRT3*M_SQRT2*R*(kx+ky)*(kx-ky)*g2+I*M_SQRT3*M_SQRT2*R*kx*ky*2.0*g3);
    setmat(h, 5, 11, 2*M_SQRT3*M_SQRT2*g3*R*kz*kp);
    setmat(h, 5, 13, -2*M_SQRT3*g3*R*kz*kp);
    
    setmat (h, 6, 4, conj(getmat(h, 4, 6)));
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    setmat (h, 6, 11, -6.0*g3*R*kz*km);
    setmat (h, 6, 12, getmat(h,5,13));
    
    setmat (h, 11, 5, conj(getmat(h, 5, 11)));
    setmat (h, 11, 6, conj(getmat(h, 6, 11)));
    setmat (h, 11, 12, M_SQRT3*R*(kx+ky)*(ky-kx)*g2-I*M_SQRT3*R*kx*ky*2.0*g3);
    setmat (h, 11, 13, getmat(h,4,6));
    
    setmat (h, 12, 4, conj(getmat(h, 4, 12)));
    setmat (h, 12, 6, conj(getmat(h, 6, 12)));
    setmat (h, 12, 11, conj(getmat(h, 11, 12)));
    setmat (h, 12, 13, M_SQRT3*M_SQRT2*R*(kx+ky)*(ky-kx)*g2+I*M_SQRT3*M_SQRT2*R*kx*ky*2.0*g3);
    
    setmat (h, 13, 4, conj(getmat(h, 4, 13)));
    setmat (h, 13, 5, conj(getmat(h, 5, 13)));
    setmat (h, 13, 11, conj(getmat(h, 11, 13)));
    setmat (h, 13, 12, conj(getmat(h, 12, 13)));

	return h;
}

static
complex double * d2Hdx2_14 () {
    complex double * h = complex_double_array_calloc (196);
    
    double EP=(material->P0)*(material->P0)/R;
	double EQ=(material->Q)*(material->Q)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg))-EQ/(3*(material->E0p))-EQ/(3*((material->E0p)+(material->D0p)));
    double g2=(material->G2L)-EP/(6*(material->Eg))+EQ/(6*(material->E0p));
    
    complex double a = 2*R;
    setmat (h, 0, 0, a);
    setmat (h, 1, 1, a);
    setmat (h, 7, 7, a);
    setmat (h, 8, 8, a);
    setmat (h, 2, 2, a);
    setmat (h, 9, 9, a);
    
    a = 2*R*2*(material->F)+2*R;  // conduction band
    setmat(h, 3, 3, a);
    setmat(h, 10, 10, a);
    
    a=-2*R*(g1-g2);
    setmat(h, 4, 4, a);
    setmat(h, 11, 11, a);
    
    a=-2*R*(g1+g2);
    setmat(h, 5, 5, a);
    setmat(h, 12, 12, a);
    
    a=-2*R*g1;
    setmat(h, 6, 6, a);
    setmat(h, 13, 13, a);
        
    setmat(h, 4, 5, M_SQRT3*2*R*g2);
    setmat(h, 4, 6, -2*R*M_SQRT2*g2);
    
    setmat(h, 5, 4, conj(getmat(h,4,5)));
    setmat(h, 5, 6, M_SQRT3*M_SQRT2*2*R*g2);
    
    setmat (h, 6, 4, conj(getmat(h, 4, 6)));
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    
    setmat (h, 11, 12, -2*M_SQRT3*R*g2);
    setmat (h, 11, 13, getmat(h,4,6));
    
    setmat (h, 12, 11, conj(getmat(h, 11, 12)));
    setmat (h, 12, 13, -2*M_SQRT3*M_SQRT2*R*g2);
    
    setmat (h, 13, 11, conj(getmat(h, 11, 13)));
    setmat (h, 13, 12, conj(getmat(h, 12, 13)));
    
    return h;
}

static
complex double * d2Hdy2_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    double EP=(material->P0)*(material->P0)/R;
	double EQ=(material->Q)*(material->Q)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg))-EQ/(3*(material->E0p))-EQ/(3*((material->E0p)+(material->D0p)));
    double g2=(material->G2L)-EP/(6*(material->Eg))+EQ/(6*(material->E0p));
    
    complex double a = 2*R;
    setmat (h, 0, 0, a);
    setmat (h, 1, 1, a);
    setmat (h, 7, 7, a);
    setmat (h, 8, 8, a);
    setmat (h, 2, 2, a);
    setmat (h, 9, 9, a);
    
    a = 2*R*2*(material->F)+2*R;  // conduction band
    setmat(h, 3, 3, a);
    setmat(h, 10, 10, a);
    
    a=-2*R*(g1-g2);
    setmat(h, 4, 4, a);
    setmat(h, 11, 11, a);
    
    a=-2*R*(g1+g2);
    setmat(h, 5, 5, a);
    setmat(h, 12, 12, a);
    
    a=-2*R*g1;
    setmat(h, 6, 6, a);
    setmat(h, 13, 13, a);
    
    setmat(h, 4, 5, -2*M_SQRT3*R*g2);
    setmat(h, 4, 6, -2*M_SQRT2*g2*R);
    
    setmat(h, 5, 4, conj(getmat(h,4,5)));
    setmat(h, 5, 6, -2*M_SQRT3*M_SQRT2*R*g2);
    
    setmat (h, 6, 4, conj(getmat(h, 4, 6)));
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    
    setmat (h, 11, 12, 2*M_SQRT3*R*g2);
    setmat (h, 11, 13, getmat(h,4,6));
    
    setmat (h, 12, 11, conj(getmat(h, 11, 12)));
    setmat (h, 12, 13, 2*M_SQRT3*M_SQRT2*R*g2);
    
    setmat (h, 13, 11, conj(getmat(h, 11, 13)));
    setmat (h, 13, 12, conj(getmat(h, 12, 13)));
    
    return h;
}

static
complex double * d2Hdz2_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    double EP=(material->P0)*(material->P0)/R;
	double EQ=(material->Q)*(material->Q)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg))-EQ/(3*(material->E0p))-EQ/(3*((material->E0p)+(material->D0p)));
    double g2=(material->G2L)-EP/(6*(material->Eg))+EQ/(6*(material->E0p));
    
    complex double a = 2*R;
    setmat (h, 0, 0, a);
    setmat (h, 1, 1, a);
    setmat (h, 7, 7, a);
    setmat (h, 8, 8, a);
    setmat (h, 2, 2, a);
    setmat (h, 9, 9, a);
    
    a = 2*R*2*(material->F)+2*R;  // conduction band
    setmat(h, 3, 3, a);
    setmat(h, 10, 10, a);
    
    a=-2*R*(g1-g2)-6*R*g2;
    setmat(h, 4, 4, a);
    setmat(h, 11, 11, a);
    
    a=-2*R*(g1+g2)+6*R*g2;
    setmat(h, 5, 5, a);
    setmat(h, 12, 12, a);
    
    a=-2*R*g1;
    setmat(h, 6, 6, a);
    setmat(h, 13, 13, a);
    
    setmat(h, 4, 6, 4*R*M_SQRT2*g2);
        
    setmat (h, 6, 4, conj(getmat(h, 4, 6)));
    
    setmat (h, 11, 13, getmat(h,4,6));
        
    setmat (h, 13, 11, conj(getmat(h, 11, 13)));
    
    return h;
}

static
complex double * d2Hdxdy_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    double EP=(material->P0)*(material->P0)/R;
    double EQ=(material->Q)*(material->Q)/R;
    double g3=(material->G3L)-EP/(6*(material->Eg))-EQ/(6*(material->E0p));
    
    setmat(h, 4, 5, -I*M_SQRT3*R*2.0*g3);
    
    setmat(h, 5, 4, conj(getmat(h,4,5)));
    setmat(h, 5, 6, I*M_SQRT3*M_SQRT2*R*2.0*g3);
    
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    
    setmat (h, 11, 12, -I*M_SQRT3*R*2.0*g3);
    
    setmat (h, 12, 11, conj(getmat(h, 11, 12)));
    setmat (h, 12, 13, I*M_SQRT3*M_SQRT2*R*2.0*g3);
    
    setmat (h, 13, 12, conj(getmat(h, 12, 13)));
    
    return h;
}

static
complex double * d2Hdxdz_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    double EP=(material->P0)*(material->P0)/R;
    double EQ=(material->Q)*(material->Q)/R;
    double g3=(material->G3L)-EP/(6*(material->Eg))-EQ/(6*(material->E0p));
    
    setmat(h, 4, 12, 2*M_SQRT3*g3*R);
    setmat(h, 4, 13, 6*g3*R/M_SQRT2);
    
    setmat(h, 5, 11, 2*M_SQRT3*g3*R);
    setmat(h, 5, 13, -2*M_SQRT3/M_SQRT2*g3*R);
    
    setmat (h, 6, 11, -6.0*g3*R/M_SQRT2);
    setmat (h, 6, 12, getmat(h, 5, 13));
    
    setmat (h, 11, 5, conj(getmat(h, 5, 11)));
    setmat (h, 11, 6, conj(getmat(h, 6, 11)));
    
    setmat (h, 12, 4, conj(getmat(h, 4, 12)));
    setmat (h, 12, 6, conj(getmat(h, 6, 12)));
    
    setmat (h, 13, 4, conj(getmat(h, 4, 13)));
    setmat (h, 13, 5, conj(getmat(h, 5, 13)));
    
    return h;
}

static
complex double * d2Hdydz_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    double EP=(material->P0)*(material->P0)/R;
    double EQ=(material->Q)*(material->Q)/R;
    double g3=(material->G3L)-EP/(6*(material->Eg))-EQ/(6*(material->E0p));
    
    setmat(h, 4, 12, 2*M_SQRT3*g3*R*I);
    setmat(h, 4, 13, 6*g3*R*(-I)/M_SQRT2);
    
    setmat(h, 5, 11, 2*M_SQRT3*g3*R*I);
    setmat(h, 5, 13, -2*M_SQRT3/M_SQRT2*g3*R*I);
    
    setmat (h, 6, 11, -6.0*g3*R*(-I)/M_SQRT2);
    setmat (h, 6, 12, getmat(h,5,13));
    
    setmat (h, 11, 5, conj(getmat(h, 5, 11)));
    setmat (h, 11, 6, conj(getmat(h, 6, 11)));
    
    setmat (h, 12, 4, conj(getmat(h, 4, 12)));
    setmat (h, 12, 6, conj(getmat(h, 6, 12)));
    
    setmat (h, 13, 4, conj(getmat(h, 4, 13)));
    setmat (h, 13, 5, conj(getmat(h, 5, 13)));
    
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
	double ky = k->x[1];
	double kz = k->x[2];
	
	px = d2Hdx2_14();
    scale_then_add_complex_double_arrays (h, px, &kx, 196);
    fftw_free(px);
    
    px = d2Hdxdy_14();
    scale_then_add_complex_double_arrays (h, px, &ky, 196);
    fftw_free(px);
    
    px = d2Hdxdz_14();
    scale_then_add_complex_double_arrays (h, px, &kz, 196);
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
    
    double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
	py = d2Hdxdy_14();
    scale_then_add_complex_double_arrays (h, py, &kx, 196);
    fftw_free(py);
    
    py = d2Hdy2_14();
    scale_then_add_complex_double_arrays (h, py, &ky, 196);
    fftw_free(py);
    
    py = d2Hdydz_14();
    scale_then_add_complex_double_arrays (h, py, &kz, 196);
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
    
    double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
	pz = d2Hdxdz_14();
    scale_then_add_complex_double_arrays (h, pz, &kx, 196);
    fftw_free(pz);
    
    pz = d2Hdydz_14();
    scale_then_add_complex_double_arrays (h, pz, &ky, 196);
    fftw_free(pz);
    
    pz = d2Hdz2_14();
    scale_then_add_complex_double_arrays (h, pz, &kz, 196);
    fftw_free(pz);
    
    return h;
}

Model * zincblende14ns_new ()
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
	
	zincblende->H = H14ns;
	
	zincblende->dHdx = dHdx14;
	zincblende->dHdy = dHdy14;
	zincblende->dHdz = dHdz14;
	
	zincblende->d2Hdx2 = d2Hdx2_14 ();
	zincblende->d2Hdy2 = d2Hdy2_14 ();
	zincblende->d2Hdz2 = d2Hdz2_14 ();
	zincblende->d2Hdxdy = d2Hdxdy_14 ();
	zincblende->d2Hdxdz = d2Hdxdz_14 ();
	zincblende->d2Hdydz = d2Hdydz_14 ();
    
    bundle14(zincblende);
	
	return zincblende;
}
