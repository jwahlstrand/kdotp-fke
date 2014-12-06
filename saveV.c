#include <complex.h>
#include <glib.h>
#include "k-dot-p.h"
#include "matrix-element-pert.h"
#include "three-vector.h"
#include "gamma-spline.h"

static gchar ** kperpn = NULL;
static gchar * prefix = NULL;
static gchar * dir = NULL;
double kcalc = 0.1;

extern Model *model;
extern Material *material;

static GOptionEntry entries[] = 
{
  { "output", 'o', 0, G_OPTION_ARG_FILENAME, &prefix, "Output prefix (default is V)", NULL },
  { "dir", 'r', 0, G_OPTION_ARG_STRING, &dir, "Integration direction (default is 0_0_1)", NULL },
  { "kpar", 'c', 0, G_OPTION_ARG_DOUBLE, &kcalc, "k_parallel to output (default is 0.1)", NULL},
  { G_OPTION_REMAINING, 'k', 0, G_OPTION_ARG_STRING_ARRAY, &kperpn, "kperp (default is 0 0 0)", NULL },
  { NULL }
};

void spew_file_double_array (char * filename, double * arr, int size1) {
    int m,n;
	FILE * f = fopen(filename, "w");
    for (m=0;m<size1;m++) {
        fprintf(f, "%1.5e\n",arr[m]);
    }
    fprintf(f,"\n");
	fclose(f);
}

void spew_file_complex_array (char * filename, complex double * arr, int size1, int size2) {
    int m,n;
	FILE * f = fopen(filename, "w");
    for (m=0;m<size1;m++) {
        for (n=0;n<size2;n++) {
            fprintf(f, "%1.5e %1.5e ",creal(arr[m+size1*n]),cimag(arr[m+size1*n]));
        }
        fprintf(f, "\n");
    }
    fprintf(f,"\n");
	fclose(f);
}

int main (int argc, char *argv[]) {
    GError *error = NULL;
    GOptionContext *context;
    
    context = g_option_context_new ("-- kperpe kperpf kperpg");
    g_option_context_add_main_entries (context, entries, NULL);
    GOptionGroup * kpgroup = k_dot_p_get_option_group();
    g_option_context_add_group (context, kpgroup);
    g_option_context_parse (context, &argc, &argv, &error);
    
    g_option_context_free(context);
    k_dot_p_init();
    
    ThreeVector *kperp =three_vector_from_string_array(kperpn);
    
    ThreeVector *kdir;
    if(dir!=NULL)
	  kdir = dir_from_string (dir);
    else {
      gchar *dstring = g_strdup("0_0_1");
      kdir = dir_from_string(dstring);
    }
    
    if(prefix==NULL)
      prefix="V";
    
    GList * list = calc_w_phi_coeffs (kperp, kdir);  // calculate list of V at kperp as a function of kpar
    GammaSpline * gs = gamma_spline_new (g_list_length(list));
    int NB = model->total_bands;
    
    double * en = gamma_spline_energy(gs,list,kcalc);
    char filename[100];
    
    sprintf(filename,"%s_en",prefix);
    spew_file_double_array(filename,en,NB);
    
    complex double * Vx = gamma_spline_V(gs,list,kcalc,0);
    
    sprintf(filename,"%s_x",prefix);
    spew_file_complex_array(filename,Vx,NB,NB);
    
    complex double * Vy = gamma_spline_V(gs,list,kcalc,1);
    
    sprintf(filename,"%s_y",prefix);
    spew_file_complex_array(filename,Vy,NB,NB);
    
    complex double * Vz = gamma_spline_V(gs,list,kcalc,2);
    sprintf(filename,"%s_z",prefix);
    spew_file_complex_array(filename,Vz,NB,NB);
    
    gamma_spline_free(gs);
    matrix_element_list_free (list);
    
    double tcp = coeffs_pert_get_elapsed_time();
    printf("coeffs_pert: %f\n",tcp);
    
    double tgs = gamma_spline_get_elapsed_time();
    printf("gamma_spline: %f\n",tgs);
    
    return 0;
}
