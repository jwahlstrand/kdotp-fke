#include "zincblende.h"
#include <gsl/gsl_math.h>

#define epsilon 0.0//1e-10
#define setmat(mat,i,j,val) mat[(j)+(i)*14]=(val)
#define getmat(mat,i,j) mat[(j)+(i)*14]

extern Material * material;
extern Model * model;

static
complex double * H14 (const ThreeVector * k)
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
    double E0 = -Eg;
    double E1 = material->E0p - material->Eg;
    
    double G0 = -Eg-material->D0;
    double G1 = E1+material->D0p;
    
    double db = material->db;
    double Ck = material->Ck;
    
    scale_then_add_complex_double_arrays (h, model->Px, &kx, 196);
    scale_then_add_complex_double_arrays (h, model->Py, &ky, 196);
    scale_then_add_complex_double_arrays (h, model->Pz, &kz, 196);
    
    double G1p=(G1+E0)/2+sqrt((G1-E0)*(G1-E0)/4-db*db/9);
    double E0pp=(G1+E0)/2-sqrt((G1-E0)*(G1-E0)/4-db*db/9);
    double G0p=(G0+E1)/2-sqrt((G0-E1)*(G0-E1)/4-4*db*db/9);
    double E1p=(G0+E1)/2+sqrt((G0-E1)*(G0-E1)/4-4*db*db/9);
    
    /*printf("E0 is %f, ",E0pp);
    printf("E1 is %f, ",E1p);
    printf("G0 is %f, ",G0p);
    printf("G1 is %f\n", G1p);*/
	
	/* diagonal elements */
	
	complex double a = G1p+Ek;
    setmat (h, 0, 0, a);
    setmat (h, 1, 1, a+epsilon);
    setmat (h, 7, 7, a);
    setmat (h, 8, 8, a+epsilon);
    
    a = E1p+Ek;
    setmat (h, 2, 2, a);
    setmat (h, 9, 9, a+epsilon);
            
    a = Ek*2*(material->F)+Ek;  // conduction band
    setmat(h, 3, 3, a);
    setmat(h, 10, 10, a+epsilon);
    
    a=E0pp-Ek*(g1-g2)-3*Ez*g2;
    setmat(h, 4, 4, a);
    setmat(h, 11, 11, a+epsilon);
    
    a=E0pp-Ek*(g1+g2)+3*Ez*g2;
    setmat(h, 5, 5, a);
    setmat(h, 12, 12, a+epsilon);
    
    a=G0p-Ek*g1;
    setmat(h, 6, 6, a);
    setmat(h, 13, 13, a+epsilon);
    
    complex double kp = (kx + I*ky)/M_SQRT2;
    complex double km=conj(kp);
    
    setmat (h, 0, 4, db/3.0);
    
    setmat (h, 1, 5, conj(getmat(h, 0, 4)));
    
    setmat (h, 2, 6, -2*db/3.0);
    
    setmat(h, 4, 0, conj(getmat(h, 0, 4)));
    setmat(h, 4, 5, M_SQRT3*R*(kx+ky)*(kx-ky)*g2-Ck*kz-I*M_SQRT3*R*kx*ky*2.0*g3);
    setmat(h, 4, 6, M_SQRT2*g2*E2zmxy);
    setmat(h, 4, 11, -M_SQRT3/M_SQRT2*Ck*kp);
    setmat(h, 4, 12, 2*M_SQRT3*M_SQRT2*g3*R*kz*kp-Ck*km/M_SQRT2);
    setmat(h, 4, 13, 6*g3*R*kz*km);
    
    setmat(h, 5, 1, conj(getmat(h,1,5)));
    setmat(h, 5, 4, conj(getmat(h,4,5)));
    setmat(h, 5, 6, M_SQRT3*M_SQRT2*R*(kx+ky)*(kx-ky)*g2+I*M_SQRT3*M_SQRT2*R*kx*ky*2.0*g3);
    setmat(h, 5, 11, 2*M_SQRT3*M_SQRT2*g3*R*kz*kp+Ck*km/M_SQRT2);
    setmat(h, 5, 12, -M_SQRT3/M_SQRT2*Ck*kp);
    setmat(h, 5, 13, -2*M_SQRT3*g3*R*kz*kp);
    
    setmat (h, 6, 2, conj(getmat(h, 2, 6)));
    setmat (h, 6, 4, conj(getmat(h, 4, 6)));
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    setmat (h, 6, 11, -6.0*g3*R*kz*km);
    setmat (h, 6, 12, getmat(h,5,13));
    
    setmat (h, 7, 11, getmat(h, 0, 4));
    
    setmat (h, 8, 12, getmat(h, 0, 4));
    
    setmat (h, 9, 13, getmat(h, 2, 6));
    
    setmat (h, 11, 4, conj(getmat(h, 4, 11)));
    setmat (h, 11, 5, conj(getmat(h, 5, 11)));
    setmat (h, 11, 6, conj(getmat(h, 6, 11)));
    setmat (h, 11, 7, conj(getmat(h, 7, 11)));
    setmat (h, 11, 12, M_SQRT3*R*(kx+ky)*(ky-kx)*g2-Ck*kz-I*M_SQRT3*R*kx*ky*2.0*g3);
    setmat (h, 11, 13, getmat(h,4,6));
    
    setmat (h, 12, 4, conj(getmat(h, 4, 12)));
    setmat (h, 12, 5, conj(getmat(h, 5, 12)));
    setmat (h, 12, 6, conj(getmat(h, 6, 12)));
    setmat (h, 12, 8, conj(getmat(h, 8, 12)));
    setmat (h, 12, 11, conj(getmat(h, 11, 12)));
    setmat (h, 12, 13, M_SQRT3*M_SQRT2*R*(kx+ky)*(ky-kx)*g2+I*M_SQRT3*M_SQRT2*R*kx*ky*2.0*g3);
    
    setmat (h, 13, 4, conj(getmat(h, 4, 13)));
    setmat (h, 13, 5, conj(getmat(h, 5, 13)));
    setmat (h, 13, 9, conj(getmat(h, 9, 13)));
    setmat (h, 13, 11, conj(getmat(h, 11, 13)));
    setmat (h, 13, 12, conj(getmat(h, 12, 13)));

	return h;
}

static
void d2Hdx2_14i (complex double *h) {
    double EP=(material->P0)*(material->P0)/R;
	double EQ=(material->Q)*(material->Q)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg))-EQ/(3*(material->E0p))-EQ/(3*((material->E0p)+(material->D0p)));
    double g2=(material->G2L)-EP/(6*(material->Eg))+EQ/(6*(material->E0p));
    
    complex_array_zero(h,196);
    
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
}

static
complex double * d2Hdx2_14 () {
    complex double * h = complex_double_array_calloc (196);
    
    d2Hdx2_14i(h);
    
    return h;
}

static
void d2Hdy2_14i (complex double *h){
    double EP=(material->P0)*(material->P0)/R;
	double EQ=(material->Q)*(material->Q)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg))-EQ/(3*(material->E0p))-EQ/(3*((material->E0p)+(material->D0p)));
    double g2=(material->G2L)-EP/(6*(material->Eg))+EQ/(6*(material->E0p));
    
    complex_array_zero(h,196);
    
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
}

static
complex double * d2Hdy2_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    d2Hdy2_14i(h);
    
    return h;
}


static
void d2Hdz2_14i (complex double *h){
    double EP=(material->P0)*(material->P0)/R;
	double EQ=(material->Q)*(material->Q)/R;
	double g1=(material->G1L)-EP/(3*(material->Eg))-EQ/(3*(material->E0p))-EQ/(3*((material->E0p)+(material->D0p)));
    double g2=(material->G2L)-EP/(6*(material->Eg))+EQ/(6*(material->E0p));
    
    complex_array_zero(h,196);
    
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
}

static
complex double * d2Hdz2_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    d2Hdz2_14i(h);
    
    return h;
}

static
void d2Hdxdy_14i (complex double *h){
    double EP=(material->P0)*(material->P0)/R;
    double EQ=(material->Q)*(material->Q)/R;
    double g3=(material->G3L)-EP/(6*(material->Eg))-EQ/(6*(material->E0p));
    
    complex_array_zero(h,196);
    
    setmat(h, 4, 5, -I*M_SQRT3*R*2.0*g3);
    
    setmat(h, 5, 4, conj(getmat(h,4,5)));
    setmat(h, 5, 6, I*M_SQRT3*M_SQRT2*R*2.0*g3);
    
    setmat (h, 6, 5, conj(getmat(h, 5, 6)));
    
    setmat (h, 11, 12, -I*M_SQRT3*R*2.0*g3);
    
    setmat (h, 12, 11, conj(getmat(h, 11, 12)));
    setmat (h, 12, 13, I*M_SQRT3*M_SQRT2*R*2.0*g3);
    
    setmat (h, 13, 12, conj(getmat(h, 12, 13)));
}

static
complex double * d2Hdxdy_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    d2Hdxdy_14i(h);
    
    return h;
}

static
void d2Hdxdz_14i (complex double *h){
    double EP=(material->P0)*(material->P0)/R;
    double EQ=(material->Q)*(material->Q)/R;
    double g3=(material->G3L)-EP/(6*(material->Eg))-EQ/(6*(material->E0p));
    
    complex_array_zero(h,196);
    
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
}

static
complex double * d2Hdxdz_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    d2Hdxdz_14i(h);
    
    return h;
}

static
void d2Hdydz_14i (complex double *h){
    double EP=(material->P0)*(material->P0)/R;
    double EQ=(material->Q)*(material->Q)/R;
    double g3=(material->G3L)-EP/(6*(material->Eg))-EQ/(6*(material->E0p));
    
    complex_array_zero(h,196);
    
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
}

static
complex double * d2Hdydz_14 (){
    complex double * h = complex_double_array_calloc (196);
    
    d2Hdydz_14i(h);
    
    return h;
}

static
void dHdx14 (complex double * h, const ThreeVector * k)
{
    complex_array_zero(h,196);
    add_complex_double_arrays (h, model->Px, 196);
    
    double Ck = material->Ck;
    
    double kx = k->x[0];
    double ky = k->x[1];
    double kz = k->x[2];
	
    setmat(h, 4, 11, -M_SQRT3*Ck/2.0);
    setmat(h, 4, 12, -Ck/2);
    
    setmat(h, 5, 11, Ck/2);
    setmat(h, 5, 12, -M_SQRT3*Ck/2.0);
    
    setmat(h, 11, 4, conj(getmat(h, 4, 11)));
    setmat(h, 11, 5, conj(getmat(h, 5, 11)));
    
    setmat(h, 12, 4, conj(getmat(h, 4, 12)));
    setmat(h, 12, 5, conj(getmat(h, 5, 12)));
	
    complex double *px = model->buffer1;
    d2Hdx2_14i(px);
    scale_then_add_complex_double_arrays (h, px, &kx, 196);
    
    d2Hdxdy_14i(px);
    scale_then_add_complex_double_arrays (h, px, &ky, 196);
    
    d2Hdxdz_14i(px);
    scale_then_add_complex_double_arrays (h, px, &kz, 196);
}

static
void dHdy14 (complex double *h, const ThreeVector * k)
{
    complex_array_zero(h,196);
    add_complex_double_arrays (h, model->Py, 196);
    
    double Ck = material->Ck;
        
    double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
	setmat(h, 4, 11, -I*M_SQRT3*Ck/2.0);
    setmat(h, 4, 12, I*Ck/2.0);
    
    setmat(h, 5, 11, -I*Ck/2.0);
    setmat(h, 5, 12, -I*M_SQRT3*Ck/2.0);
    
    setmat(h, 11, 4, conj(getmat(h, 4, 11)));
    setmat(h, 11, 5, conj(getmat(h, 5, 11)));
    
    setmat(h, 12, 4, conj(getmat(h, 4, 12)));
    setmat(h, 12, 5, conj(getmat(h, 5, 12)));
	
	complex double *py = model->buffer1;
    d2Hdxdy_14i(py);
    scale_then_add_complex_double_arrays (h, py, &kx, 196);
    
    d2Hdy2_14i(py);
    scale_then_add_complex_double_arrays (h, py, &ky, 196);
    
    d2Hdydz_14i(py);
    scale_then_add_complex_double_arrays (h, py, &kz, 196);
}

static
void dHdz14 (complex double *h, const ThreeVector * k)
{
    complex_array_zero(h,196);
    add_complex_double_arrays (h, model->Pz, 196);
    
    double Ck = material->Ck;
        
    double kx = k->x[0];
	double ky = k->x[1];
	double kz = k->x[2];
	
	setmat (h, 4, 5, -Ck);
    setmat (h, 5, 4, -Ck);
    
    setmat (h, 11, 12, -Ck);
    setmat (h, 12, 11, -Ck);
	
	complex double *pz = model->buffer1;
    d2Hdxdz_14i(pz);
    scale_then_add_complex_double_arrays (h, pz, &kx, 196);
    
    d2Hdydz_14i(pz);
    scale_then_add_complex_double_arrays (h, pz, &ky, 196);
    
    d2Hdz2_14i(pz);
    scale_then_add_complex_double_arrays (h, pz, &kz, 196);
}

void bundle14 (Model * zincblende)
{
    	int n;
	GList * b = NULL;
	GArray * array = g_array_sized_new (FALSE, TRUE, sizeof(int), 6);
	n=0;
	g_array_append_val (array, n);
	n=1;
	g_array_append_val (array, n);
	n=2;
	g_array_append_val (array, n);
	n=3;
	g_array_append_val (array, n);
	n=4;
	g_array_append_val (array, n);
	n=5;
	g_array_append_val (array, n);
	
	Bundle * b1 = g_slice_new(Bundle);
	b1->list = array;
	b1->first = 0;
	
	b=g_list_append (b, b1);
	
	GArray * array2 = g_array_sized_new (FALSE, TRUE, sizeof(int), 2);
	n=6;
	g_array_append_val (array2, n);
	n=7;
	g_array_append_val (array2, n);
	
	Bundle * b2 = g_slice_new(Bundle);
	b2->list = array2;
	b2->first = 6;
	
	b=g_list_append (b, b2);

#if UCBANDS
    GArray * array3 = g_array_sized_new (FALSE, TRUE, sizeof(int), 6);
	n=8;
	g_array_append_val (array3, n);
	n=9;
	g_array_append_val (array3, n);
	n=10;
	g_array_append_val (array3, n);
	n=11;
	g_array_append_val (array3, n);
	n=12;
	g_array_append_val (array3, n);
	n=13;
	g_array_append_val (array3, n);
	
	Bundle * b3 = g_slice_new(Bundle);
	b3->list = array3;
	b3->first = 8;
	
	b=g_list_append (b, b3);
#endif	
	zincblende->bundles = b;

}

Model * zincblende14_new ()
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
    
    zincblende->buffer1 = complex_double_array_calloc (196);
	
	zincblende->H = H14;
	
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

