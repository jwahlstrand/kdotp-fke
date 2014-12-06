#include <string.h>
#include "k-dot-p.h"
#include "zincblende.h"

#define OPTIMIZE 1
#define SUPEROPTIMIZE 1

static gchar * modelname = "zincblende14";
static gchar * materialname = "GaAs";

extern Material * material = NULL;
extern Model * model = NULL;

#define HADER 0 /* use Hader parameters to compare to the old 8-band FKE paper */

static GOptionEntry kdotp_entries[] = 
{
  { "model", 'm', 0, G_OPTION_ARG_STRING, &modelname, "Band structure model to use (default is zincblende14) - options are zincblende2pba, zincblende14ns, zincblende14nr, and zincblende14", NULL },
  { "semiconductor", 's', 0, G_OPTION_ARG_STRING, &materialname, "Semiconductor to model (default is GaAs) - options are GaAs, GaSb, InP, InSb, and ZnSe", NULL },
  { NULL }
};

complex double  * complex_double_array_calloc (size_t n) {
    complex double  * array = (complex double *) fftw_malloc(sizeof(complex double) * n);
    memset(array, 0, sizeof(complex double )*n);
    return array;
}

void add_complex_double_arrays (complex double * darray, complex double * sarray, size_t n)
{
    size_t i;
    for (i=0;i<n;i++) {
        darray[i]+=sarray[i];
    }
}

void scale_then_add_complex_double_arrays (complex double * darray, complex double * sarray, double * factor, size_t n)
{
    size_t i;
    for (i=0;i<n;i++) {
        darray[i]+=sarray[i]*(*factor);
    }
}

/* B = U^H A U for NxN square matrices */

void matrix_transform (complex double * B, const complex double * U, const complex double * A, size_t N)
{
    size_t m, n, q, p;
    
    // B should be zeroed!
    
    /* assume that B and A are both Hermitian */
#if SUPEROPTIMIZE
    complex double litlsum1,litlsum2;
    for (n=0;n<N;n++) {
        for (q=0;q<N;q+=2) {
            litlsum1 = 0.0;
            litlsum2 = 0.0;
            for (p=0;p<N;p++) {
                litlsum1+=U[p*N+n]*A[q*N+p];
                litlsum2+=U[p*N+n]*A[(q+1)*N+p];
            }
            for (m=0;m<=n;m++) {
                B[m*N+n]+=conj(U[q*N+m])*litlsum1;
                B[m*N+n]+=conj(U[(q+1)*N+m])*litlsum2;
            }
        }
    }
    for (n=0;n<N;n++) {
        for (m=n+1;m<N;m++) {
            B[m*N+n]=conj(B[n*N+m]);
        }
    }
#elif OPTIMIZE
    complex double litlsum1, litlsum2;

    for (m=0;m<N;m+=2) {
        for (n=0;n<N;n++) {
            if (n>=m) {
                for (q=0;q<N;q+=2) {
                    litlsum1 = 0.0;
                    litlsum2 = 0.0;
                    for (p=0;p<N;p++) {
                        litlsum1+=U[p*N+n]*A[q*N+p];
                        litlsum2+=U[p*N+n]*A[(q+1)*N+p];
                    }
                    B[m*N+n]+=conj(U[q*N+m])*litlsum1;
                    B[m*N+n]+=conj(U[(q+1)*N+m])*litlsum2;
                    B[(m+1)*N+n]+=conj(U[q*N+(m+1)])*litlsum1;
                    B[(m+1)*N+n]+=conj(U[(q+1)*N+(m+1)])*litlsum2;
                }
            }
            else {
                B[m*N+n]=conj(B[n*N+m]);
                B[(m+1)*N+n]=conj(B[n*N+(m+1)]);
            }
        }
    }
#else
    complex double litlsum;

    for (m=0;m<N;m++) {
        for (n=0;n<N;n++) {
            if (n>=m) {
                for (q=0;q<N;q++) {
                    litlsum = 0.0;
                    for (p=0;p<N;p++) {
                        litlsum+=U[p*N+n]*A[q*N+p];
                    }
                    B[m*N+n]+=conj(U[q*N+m])*litlsum;
                }
            }
            else {
                B[m*N+n]=conj(B[n*N+m]);
            }
        }
    }
#endif
}

void matrix_transform3 (complex double * Bx, complex double * By, complex double * Bz, const complex double * U, const complex double * Ax, const complex double * Ay, const complex double * Az, size_t N)
{
    size_t m, n, q, p;
    complex double Cconj, litlsumx,litlsumy,litlsumz;
    
    /* assume that B and A are both Hermitian */
#if SUPEROPTIMIZE
    for (n=0;n<N;n++) {
        for (q=0;q<N;q++) {
            litlsumx = 0.0;
            litlsumy = 0.0;
            litlsumz = 0.0;
            for (p=0;p<N;p++) {
                litlsumx+=U[p*N+n]*Ax[q*N+p];
                litlsumy+=U[p*N+n]*Ay[q*N+p];
                litlsumz+=U[p*N+n]*Az[q*N+p];
            }
            for (m=0;m<=n;m++) {
                Cconj = conj(U[q*N+m]);
                Bx[m*N+n]+=Cconj*litlsumx;
                By[m*N+n]+=Cconj*litlsumy;
                Bz[m*N+n]+=Cconj*litlsumz;
            }
        }
    }
    for (n=0;n<N;n++) {
        for (m=n+1;m<N;m++) {
            Bx[m*N+n]=conj(Bx[n*N+m]);
            By[m*N+n]=conj(By[n*N+m]);
            Bz[m*N+n]=conj(Bz[n*N+m]);
        }
    }
#elif OPTIMIZE
    for (m=0;m<N;m++) {
        for (n=0;n<N;n++) {
            if (n>=m) {
                for (q=0;q<N;q++) {
                    Cconj = conj(U[q*N+m]);
                    litlsumx = 0.0;
                    litlsumy = 0.0;
                    litlsumz = 0.0;
                    for (p=0;p<N;p++) {
                        litlsumx+=U[p*N+n]*Ax[q*N+p];
                        litlsumy+=U[p*N+n]*Ay[q*N+p];
                        litlsumz+=U[p*N+n]*Az[q*N+p];
                    }
                    Bx[m*N+n]+=Cconj*litlsumx;
                    By[m*N+n]+=Cconj*litlsumy;
                    Bz[m*N+n]+=Cconj*litlsumz;
                }
            }
            else {
                Bx[m*N+n]=conj(Bx[n*N+m]);
                By[m*N+n]=conj(By[n*N+m]);
                Bz[m*N+n]=conj(Bz[n*N+m]);
            }
        }
    }
#else
    for (m=0;m<N;m++) {
        for (n=m;n<N;n++) {
            for (q=0;q<N;q++) {
                Cconj = conj(U[q*N+m]);
                litlsumx = 0.0;
                litlsumy = 0.0;
                litlsumz = 0.0;
                for (p=0;p<N;p++) {
                    litlsumx+=U[p*N+n]*Ax[q*N+p];
                    litlsumy+=U[p*N+n]*Ay[q*N+p];
                    litlsumz+=U[p*N+n]*Az[q*N+p];
                }
                Bx[m*N+n]+=Cconj*litlsumx;
                By[m*N+n]+=Cconj*litlsumy;
                Bz[m*N+n]+=Cconj*litlsumz;
            }
        }
        for(n=0;n<m;n++) {
            Bx[m*N+n]=conj(Bx[n*N+m]);
            By[m*N+n]=conj(By[n*N+m]);
            Bz[m*N+n]=conj(Bz[n*N+m]);
        }
    }
#endif
}

void k_dot_p_init () {
	/* set material, model from materialname, modelname */
	
	material = material_new_by_name (materialname);
    model = model_new_by_name (modelname);
    
    
}

void k_dot_p_cleanup () {
    g_slice_free (Material, material);
    model_free(model);
    fftw_cleanup();
}

GOptionGroup * k_dot_p_get_option_group () {
    GOptionGroup * gog = g_option_group_new ("kp","General options:", "Show k.p options", NULL, NULL);
    g_option_group_add_entries (gog, kdotp_entries);
    return gog;
}

Material * material_new_by_name (gchar * name) {
    if(!strcmp(name,"GaAs")) {
    	Material * GaAs = g_slice_new (Material);
    	GaAs->P0 = 10.30;
    	GaAs->Q = 7.70;
    	GaAs->P0p = 3.00;
    	GaAs->Eg = 1.519;
    	GaAs->E0p = 4.488;
    	GaAs->D0 = 0.341;
    	GaAs->D0p = 0.171;
    	GaAs->G1L = 7.797;
    	GaAs->G2L = 2.458;
    	GaAs->G3L = 3.299;
    	GaAs->F = -1.055;
    	GaAs->db = -0.061;
    	GaAs->Ck = -0.0034;
    	GaAs->D1 = 0.081;
    	GaAs->D2 = 0.0;
    	GaAs->D3 = 0.211;
    	GaAs->E6v = 12.55;
    	GaAs->E6u = 8.56;
    	GaAs->E3d = 10.17;
    	GaAs->E5d = 11.89;
    	GaAs->E6q = 13.64;
    	GaAs->P0_30 = 9.232;
    	GaAs->Qc = 7.998;
    	GaAs->P3 = 4.328;
    	GaAs->P2 = 4.891;
    	GaAs->Ps = 3.045;
    	GaAs->Pd = 0.195;
    	GaAs->Qcd = 4.068;
    	GaAs->P3d = 5.819;
    	GaAs->P2d = 9.392;
    	GaAs->Pu = 8.648;
    	GaAs->Pm1 = 0.5;
    #if HADER
        GaAs->P0 = 9.98;
        GaAs->G1L = 6.85;
        GaAs->G2L = 2.1;
        GaAs->G3L = 2.90;
    #endif
    	return GaAs;
    }
    if(!strcmp(name,"ZnSe")) {
        Material * ZnSe = g_slice_new (Material);
        ZnSe->P0 = 10.628;
        ZnSe->Q = 9.845;
        ZnSe->P0p = 9.165;
        ZnSe->Eg = 2.82;
        ZnSe->E0p = 7.33;
        ZnSe->D0 = 0.403;
        ZnSe->D0p = 0.09;
        ZnSe->G1L = 4.3;
        ZnSe->G2L = 1.14;
        ZnSe->G3L = 1.84;
        ZnSe->F = 0.0;
        ZnSe->db = -0.238;
        ZnSe->Ck = -0.014;
        return ZnSe;
    }
    if(!strcmp(name,"InP")) {
        Material * InP = g_slice_new (Material);
        InP->P0 = 8.65;
        InP->Q = 7.42;
        InP->P0p = 4.3;
        InP->Eg = 1.424;
        InP->E0p = 4.6;
        InP->D0 = 0.108;
        InP->D0p = 0.50;
        InP->G1L = 5.05;
        InP->G2L = 1.6;
        InP->G3L = 1.73;
        InP->F = 0.0;
        InP->db = 0.22;
        InP->Ck = -0.014;
        return InP;
    }
    if(!strcmp(name,"GaSb")) {
        Material * GaSb = g_slice_new (Material);
        GaSb->P0 = 9.5;
        GaSb->Q = 8.12;
        GaSb->P0p = 3.33;
        GaSb->Eg = 0.813;
        GaSb->E0p = 3.3;
        GaSb->D0 = 0.75;
        GaSb->D0p = 0.33;
        GaSb->G1L = 13.2;
        GaSb->G2L = 4.4;
        GaSb->G3L = 5.7;
        GaSb->F = 0.0;
        GaSb->db = -0.28;
        GaSb->Ck = 0.00043;
        return GaSb;
    }
    if(!strcmp(name,"InSb")) {
        Material * InSb = g_slice_new (Material);
        InSb->P0 = 9.51;
        InSb->Q = 8.22;
        InSb->P0p = 3.17;
        InSb->Eg = 0.235;
        InSb->E0p = 3.39;
        InSb->D0 = 0.803;
        InSb->D0p = 0.39;
        InSb->G1L = 40.1;
        InSb->G2L = 18.1;
        InSb->G3L = 19.2;
        InSb->F = 0.0;
        InSb->db = -0.244;
        InSb->Ck = -0.0092;
        return InSb;
    }
    g_assert_not_reached();
    return NULL;
}

Model * model_new_by_name (gchar * name) {
    g_assert (name);
    if (!strcmp(name,"zincblende2pba"))
        return zincblende2pba_new ();
    if (!strcmp(name,"zincblende14ns"))
        return zincblende14ns_new ();
    if (!strcmp(name,"zincblende14nr"))
        return zincblende14nr_new ();
    if (!strcmp(name,"zincblende14"))
        return zincblende14_new ();
    g_assert_not_reached();
    return NULL;
}

gboolean model_has_remote_params () {
    return (model->d2Hdx2 != NULL);
}

void model_free (Model * model) {
    m_free(model->Px);
    m_free(model->Py);
    m_free(model->Pz);
    if (model->d2Hdx2)
        m_free(model->d2Hdx2);
    if (model->d2Hdy2)
        m_free(model->d2Hdy2);
    if (model->d2Hdz2)
        m_free(model->d2Hdz2);
    if (model->d2Hdxdy)
        m_free(model->d2Hdxdy);
    if (model->d2Hdxdz)
        m_free(model->d2Hdxdz);
    if (model->d2Hdydz)
        m_free(model->d2Hdydz);
    g_slice_free (Model, model);
}
