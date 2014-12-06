#include <glib.h>
#include "k-dot-p.h"
#include "matrix-element-pert.h"
#include "three-vector.h"
#include "gamma-spline.h"
#include "calc-gamma.h"
#include "zener.h"

static gchar ** kperpn = NULL;
static gchar * prefix = NULL;
static gchar * dir = NULL;
double kcalc = 0.1;
double epsilon = 1.0;
int b = 0;

extern Model *model;
extern Material *material;

static GOptionEntry entries[] = 
{
  { "output", 'o', 0, G_OPTION_ARG_FILENAME, &prefix, "Output prefix (default is m)", NULL },
  { "dir", 'r', 0, G_OPTION_ARG_STRING, &dir, "Integration direction (default is 0_0_1)", NULL },
  { "kpar", 'c', 0, G_OPTION_ARG_DOUBLE, &kcalc, "k_parallel to output (default is 0.1)", NULL},
  { "band", 'b', 0, G_OPTION_ARG_INT, &b, "band of interest (determines which bundle will be used, default is 0, which is the valence band bundle)", NULL},
  { "epsilon", 'e', 0, G_OPTION_ARG_DOUBLE, &epsilon, "normalized electric field (default is 1.0)", NULL},
  { G_OPTION_REMAINING, 'k', 0, G_OPTION_ARG_STRING_ARRAY, &kperpn, "kperp (default is 0 0 0)", NULL },
  { NULL }
};

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
    
    if(prefix==NULL)
      prefix=g_strdup("m");
        
    ThreeVector *kperp = three_vector_from_string_array(kperpn);
    
    if(dir==NULL)
      dir=g_strdup("0_0_1");
    
	ThreeVector *kdir = dir_from_string (dir);
    
    GList * list = calc_w_phi_coeffs (kperp, kdir);  // calculate list of V at kperp as a function of kpar
    int NB = model->total_bands;
    
    Bundle * bundle = find_bundle (b);
    const size_t NSB = bundle->list->len;
    complex double * Lsb = zener_renormalize_bundle (list, epsilon, kdir, bundle);
    size_t Nkc = Nkc_from_epsilon(epsilon);
    double dkc = KCMAX/Nkc;
    
    size_t index = kcalc/dkc;
    
    char filename[100];
    sprintf(filename,"%s_m",prefix);
    
    spew_file_complex_array(filename, &Lsb[NSB*NSB*index],NSB,NSB);
    
    double tcp = coeffs_pert_get_elapsed_time();
    printf("coeffs_pert: %f\n",tcp);
    
    double tzn = zener_get_elapsed_time();
    printf("zener: %f\n",tzn);
    
    return 0;
}
