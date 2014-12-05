#include <glib.h>
#include <gsl/gsl_math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include "three-vector.h"

ThreeVector* three_vector_new (double xp, double y, double z)
{
    ThreeVector *tv = g_slice_new(ThreeVector);
    tv->x[0]=xp;
    tv->x[1]=y;
    tv->x[2]=z;
    return tv;
}

void three_vector_free(ThreeVector *tv) {
    g_slice_free(ThreeVector, tv);
}


double __attribute__((pure)) three_vector_length (const ThreeVector *tv)
{
    return gsl_hypot3(tv->x[0],tv->x[1],tv->x[2]);
}

ThreeVector * dir_from_string (char * string)
{
    char * p = strtok(string,"_");
    ThreeVector *k =g_slice_new0(ThreeVector);
    k->x[0]=strtod(p,NULL);
    p = strtok(NULL,"_");
    k->x[1]=strtod(p,NULL);
    p = strtok(NULL,"_}");
    k->x[2]=strtod(p,NULL);
    ThreeVector *uk = three_vector_unit(k);
    return uk;
}

ThreeVector * three_vector_unit (const ThreeVector *tv)
{
    ThreeVector *nk = g_slice_new0(ThreeVector);
	double kabs = three_vector_length (tv);
	nk->x[0]=(tv->x[0])/kabs;
	nk->x[1]=(tv->x[1])/kabs;
	nk->x[2]=(tv->x[2])/kabs;
    return nk;
}

/*******************************/

ThreeVector * three_vector_copy (const ThreeVector * k)
{
    ThreeVector * nk = g_slice_new(ThreeVector);  // free me later
    nk->x[0]=k->x[0];
    nk->x[1]=k->x[1];
    nk->x[2]=k->x[2];
    
    return nk;
}

ThreeVector * three_vector_add (const ThreeVector * k1, const ThreeVector * k2)
{
    ThreeVector * nk = g_slice_new(ThreeVector);  // free me later
    nk->x[0]=(k1->x[0]) + (k2->x[0]);
    nk->x[1]=(k1->x[1]) + (k2->x[1]);
    nk->x[2]=(k1->x[2]) + (k2->x[2]);
    
    return nk;
}

ThreeVector * three_vector_diff (const ThreeVector * k1, const ThreeVector * k2)
{
    ThreeVector * nk = g_slice_new(ThreeVector);  // free me later
    nk->x[0]=(k1->x[0]) - (k2->x[0]);
    nk->x[1]=(k1->x[1]) - (k2->x[1]);
    nk->x[2]=(k1->x[2]) - (k2->x[2]);
    
    return nk;
}

void three_vector_incr (ThreeVector *k, const ThreeVector *dk)
{
    k->x[0]+= (dk->x[0]);
    k->x[1]+= (dk->x[1]);
    k->x[2]+= (dk->x[2]);
}

void three_vector_scale (ThreeVector *k, double d)
{
    k->x[0]*= d;
    k->x[1]*= d;
    k->x[2]*= d;
}

__attribute__((pure))
double three_vector_dot (const ThreeVector *k1, const ThreeVector *k2)
{
    return (k1->x[0])*(k2->x[0])+(k1->x[1])*(k2->x[1])+(k1->x[2])*(k2->x[2]);
}

ThreeVector * three_vector_cross (const ThreeVector *k1, const ThreeVector *k2)
{
    ThreeVector * nk = g_slice_new(ThreeVector);  // free me later
    nk->x[0] = k1->x[1]*k2->x[2] - k1->x[2]*k2->x[1];
	nk->x[1] = k1->x[2]*k2->x[0] - k1->x[0]*k2->x[2];
	nk->x[2] = k1->x[0]*k2->x[1] - k1->x[1]*k2->x[0];
    return nk;
}

__attribute__((pure))
int three_vector_equal (const ThreeVector *k1, const ThreeVector *k2) {
    return (fabs(k1->x[0]-k2->x[0])<1e-6)&&(fabs(k1->x[1]-k2->x[1])<1e-6)&&(fabs(k1->x[2]-k2->x[2])<1e-6);
}

/***************/

ThreeVector * three_vector_from_string (char * text)
{
    ThreeVector *tv = g_slice_new(ThreeVector);
    double *x = tv->x;
    if (text) {
        char * p = strtok(text,",")+1;
        x[0]=strtod(p,NULL);
        p = strtok(NULL,",");
        x[1]=strtod(p,NULL);
        p = strtok(NULL,",}");
        x[2]=strtod(p,NULL);
    }
    else {
        x[0]=0.0;
        x[1]=0.0;
        x[2]=0.0;
    }
    return tv;
}

ThreeVector * three_vector_from_string_array (const char ** toks)
{
   ThreeVector *tv = g_slice_new0(ThreeVector);
   if (toks) {
        double kc = strtod(toks[0], NULL);
        if (fabs(kc)<1e-10) kc = 0;
        tv->x[0]=kc;
        kc = strtod(toks[1], NULL);
        if (fabs(kc)<1e-10) kc = 0;
        tv->x[1]=kc;
        if(toks[2]) {
            kc = strtod(toks[2], NULL);
            if (fabs(kc)<1e-10) kc = 0;
            tv->x[2]=kc;
        }
        else tv->x[2]=0.0;
    }
  return tv;
}

/*********************/
/* functions for getting basis vectors */

ThreeVector * three_vector_e (const ThreeVector *g)
{
    // this doesn't work if g==y
    if ((g->x[0]==0 && g->x[1]==1 && g->x[2]==0)) {
        return three_vector_new (0.0,0.0,1.0);
    }
	ThreeVector *y = three_vector_new (0.0,1.0,0.0);
	ThreeVector *e = three_vector_cross (y, g);
	ThreeVector *be = three_vector_unit(e);
	return be;
}

ThreeVector * three_vector_f (const ThreeVector *g, const ThreeVector *e )
{
    ThreeVector * c = three_vector_cross (g,e);
    return c;
}

double distance_between_lines (const ThreeVector *x0, const ThreeVector *x1, const ThreeVector *n0, const ThreeVector *n1) {
    double d;
    ThreeVector * x1mx0 = three_vector_diff (x1,x0);
    if (three_vector_equal (n0, n1)) {
        d = three_vector_length(x1mx0);
        three_vector_free(x1mx0);
        return d;
    }
    ThreeVector * n1cn0 = three_vector_cross (n1, n0);
    d = fabs(three_vector_dot(x1mx0,n1cn0))/fabs(three_vector_length(n1cn0));
    three_vector_free(x1mx0);
    three_vector_free(n1cn0);
    return d;
}

int trajectory_intersects_bad (const ThreeVector * kperp, const ThreeVector * g) {
    // transform kperp into x, y, z coords
    
    int intersects = FALSE;
    
    ThreeVector * origin = three_vector_new (0,0,0);
    
    //printf("kperp2 is (%e %e %e)\n", kperp->x[0],kperp->x[1],kperp->x[2]);
    
    ThreeVector * bad = three_vector_new (1,0,0);
    double d = distance_between_lines (kperp, origin, g, bad);
    //printf("distance to (100) is %e, ", d);
    three_vector_free(bad);
    
    if (d<1e-4)
        intersects = TRUE;
    
    bad = three_vector_new (0,1,0);
    d = distance_between_lines (kperp, origin, g, bad);
    //printf("distance to (010) is %e, ", d);
    three_vector_free(bad);
    
    if (d<1e-4)
        intersects = TRUE;
    
    bad = three_vector_new (0,0,1);
    d = distance_between_lines (kperp, origin, g, bad);
    //printf("distance to (001) is %e, ", d);
    three_vector_free(bad);
    
    if (d<1e-4)
        intersects = TRUE;
    
    bad = three_vector_new (1/M_SQRT3,1/M_SQRT3,1/M_SQRT3);
    d = distance_between_lines (kperp, origin, g, bad);
    //printf("distance to (111) is %e, ", d);
    three_vector_free(bad);
    
    if (d<1e-4)
        intersects = TRUE;
    
    bad = three_vector_new (-1/M_SQRT3,1/M_SQRT3,1/M_SQRT3);
    d = distance_between_lines (kperp, origin, g, bad);
    //printf("distance to (-111) is %e, ", d);
    three_vector_free(bad);
    
    if (d<1e-4)
        intersects = TRUE;
    
    bad = three_vector_new (1/M_SQRT3,-1/M_SQRT3,1/M_SQRT3);
    d = distance_between_lines (kperp, origin, g, bad);
    //printf("distance to (1-11) is %e, ", d);
    three_vector_free(bad);
    
    if (d<1e-4)
        intersects = TRUE;
    
    bad = three_vector_new (1/M_SQRT3,1/M_SQRT3,-1/M_SQRT3);
    d = distance_between_lines (kperp, origin, g, bad);
    //printf("distance to (11-1) is %e\n", d);
    three_vector_free(bad);
    
    if (d<1e-4)
        intersects = TRUE;
    
    three_vector_free(origin);
    
    return intersects;
}

/********************/
/* misc. */

void three_vector_array_project_inplace (double *P, const double * x, const double * y, const double * z, const ThreeVector * dir, size_t N)
{
    size_t m;
    memset(P, 0, sizeof(double)*N);
    if(fabs(dir->x[0])>1e-5*three_vector_length(dir)) {
      for (m=0;m<N;m++) {
        P[m] += dir->x[0]*x[m];
      }
    }
    if(fabs(dir->x[1])>1e-5*three_vector_length(dir)) {
      for (m=0;m<N;m++) {
        P[m] += dir->x[1]*y[m];
      }
    }
    if(fabs(dir->x[2])>1e-5*three_vector_length(dir)) {
      for (m=0;m<N;m++) {
        P[m] += dir->x[2]*z[m];
      }
    }
}

double * three_vector_array_project (const double * x, const double * y, const double * z, const ThreeVector * dir, size_t N)
{
    double * P = (double*) fftw_malloc(sizeof(double)*N);
    three_vector_array_project_inplace(P,x,y,z,dir,N);
    return P;
}
