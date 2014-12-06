#include <glib.h>
#include "k-dot-p.h"
#include "matrix-element-pert.h"
#include "three-vector.h"
#include "absorption.h"
#include "gamma-spline.h"
#include "box.h"
#include "integration.h"
#include "zener.h"

#define SIMPLE 1

/* twophoton */
/* calculates one-photon, QIC, and two-photon Franz-Keldysh spectra     */

static gchar ** kperpn = NULL;
static gchar * prefix = NULL; //g_strdup("tp_out");
static gchar * dir = NULL; //g_strdup("0_0_1");
static double width = 0.01;
static double maxeps = 0.5;
static int numeps = 6;
static gboolean no_two_pho = FALSE;
static gboolean dont_incl_zero = FALSE;

/* integration */
static int maxdepth = 4;

/* spectrum parameters (in eV) */
//static double min_energy = 1.0;
//static double max_energy = 2.0;

extern Model *model;
extern Material *material;

static GOptionEntry entries[] = 
{
  { "width", 'w', 0, G_OPTION_ARG_DOUBLE, &width, "Width of box to integrate (default is 0.01)", NULL},
  { "maxepsilon", 'e', 0, G_OPTION_ARG_DOUBLE, &maxeps, "Maximum epsilon to calculate (scaled value of electric field, default is 0.5", NULL},
  { "numepsilon", 'n', 0, G_OPTION_ARG_INT, &numeps, "Number of epsilons to calculate (default is 6)", NULL},
  { "zero", 'z', 0, G_OPTION_ARG_NONE, &dont_incl_zero, "Do not calculate zero-field spectrum (default is to calculate it)", NULL},
  { "maxdepth", 'd', 0, G_OPTION_ARG_INT, &maxdepth, "Maximum integration depth", NULL},
  { "output", 'o', 0, G_OPTION_ARG_FILENAME, &prefix, "Output filename (default is tp_out)", NULL },
  { "dir", 'r', 0, G_OPTION_ARG_STRING, &dir, "Integration direction (default is 0_0_1)", NULL },
  { "notwophoton", 't', 0, G_OPTION_ARG_NONE, &no_two_pho, "Do not calculate QUIC and two-photon absorption (default is no)", NULL},
  { G_OPTION_REMAINING, 'k', 0, G_OPTION_ARG_STRING_ARRAY, &kperpn, "Center of box to integrate (default is 0 0 0)", NULL },
  { NULL }
};

int main (int argc, char *argv[]) {
    GError *error = NULL;
    GOptionContext *context;
    
    GTimer * global_timer = g_timer_new();
    g_timer_start(global_timer);

    context = g_option_context_new ("-- kperpe kperpf kperpg");
    g_option_context_add_main_entries (context, entries, NULL);
    GOptionGroup * kpgroup = k_dot_p_get_option_group();
    g_option_context_add_group (context, kpgroup);
    g_option_context_parse (context, &argc, &argv, &error);
    
    g_option_context_free(context);
    k_dot_p_init();
    
    if (maxdepth<0) {
        printf("maximum depth should be 0 or greater\n");
        exit(1);
    }
    
    ThreeVector *kperp = three_vector_from_string_array(kperpn);
    
    if(dir==NULL)
      dir=g_strdup("0_0_1");
	ThreeVector *kdir = dir_from_string (dir);
	
    Box * b = new_box (kperp, kdir, width);
    
#if CAL
    g_message("exporting calibration files");
#endif
    
    OnePhoton ** array;
    
    simple_init();
    array = solve_box_simpler (NULL, b, 1e-4, maxdepth, maxeps, numeps, !dont_incl_zero);

    g_assert (array);
    
    if(prefix==NULL)
      prefix="tp_out";
	FILE * f = fopen(prefix, "wb");

    int p;
    OnePhoton * op = NULL;
    for (p=0;p<numeps;p++) {
        op = array[p];
        g_assert (op);
        one_photon_write(op,f);
	}
	fclose(f);
	
	//delete array;
	
	g_slice_free (Box, b);
    
    g_timer_stop(global_timer);
    double total_time = g_timer_elapsed(global_timer, NULL);
    double coeffs_pert_frac = coeffs_pert_get_elapsed_time()/total_time*100.0;
    double zener_frac = zener_get_elapsed_time()/total_time*100.0;
    double calc_gamma_frac = calc_gamma_get_elapsed_time()/total_time*100.0;
    double gamma_spline_frac = gamma_spline_get_elapsed_time()/total_time*100.0;
    g_message("coeffs_pert: %f%%", coeffs_pert_frac);
    g_message("calc_gamma: %f%%", calc_gamma_frac);
    g_message("zener: %f%%", zener_frac);
    g_message("gamma_spline: %f%%", gamma_spline_frac);
    g_message("other: %f%%",100.0-coeffs_pert_frac-calc_gamma_frac-zener_frac-gamma_spline_frac);
    g_message("total time elapsed: %f s",total_time);
    
    fftw_export_wisdom_to_filename("kdotp_fftw_wisdom");
    
    //k_dot_p_cleanup();
    return 0;
}
