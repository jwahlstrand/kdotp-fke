#include <glib.h>
#include "k-dot-p.h"
#include "three-vector.h"
#include "diagonalize.h"

/* savebands */
/* diagonalizes the Hamiltonian and outputs the bands along a k trajectory to a text file */

static gchar * begin_k = NULL;
static gchar * end_k = NULL;
static int steps = 100;
static gchar * filename = NULL;

extern Model *model;
extern Material *material;

static GOptionEntry entries[] = 
{
    { "begin", 'b', 0, G_OPTION_ARG_STRING, &begin_k, "Starting value in k space (default is \"{0,0,0}\")", NULL },
    { "end", 'e', 0, G_OPTION_ARG_STRING, &end_k, "Ending value in k space (default is \"{1,0,0}\")", NULL },
    { "steps", 'n', 0, G_OPTION_ARG_INT, &steps, "Steps (default is 100)", NULL},
    { "output", 'o', 0, G_OPTION_ARG_FILENAME, &filename, "Output filename (default writes to stdout)", NULL },
    { NULL }
};

int main (int argc, char *argv[]) {
    GError *error = NULL;
    GOptionContext *context;

    context = g_option_context_new ("- output band energies along a line in k space to an ASCII file");
    g_option_context_add_main_entries (context, entries, NULL);
    g_option_context_add_group (context, k_dot_p_get_option_group());
    g_option_context_parse (context, &argc, &argv, &error);
    
    g_option_context_free(context);
    
    k_dot_p_init();
    
    if (begin_k==NULL)
        begin_k = g_strdup("{0,0,0}");
    if (end_k==NULL)
        end_k = g_strdup("{1,0,0}");
    
    ThreeVector *start = three_vector_from_string(begin_k);
    ThreeVector *end = three_vector_from_string(end_k);
    
    ThreeVector *kdiff = three_vector_diff(end,start);
    three_vector_scale(kdiff,1/((double)steps));
    ThreeVector *kval = start;
    
    FILE * outstream;
    if (filename==NULL)
        outstream = stdout;
    else
        outstream = fopen (filename, "w");
    
    size_t n;
    int q;
    for (q=0;q<steps;q++) {
        complex double * H = model->H(kval);
        double * energies = eigenvalues_herm (H, model->total_bands);
        fprintf(outstream,"%e %e %e ",kval->x[0],kval->x[1],kval->x[2]);
        for (n=0;n<model->total_bands;n++) {
            fprintf(outstream, "%1.12e ", energies[n]);
        }
        fprintf(outstream,"\n");
        m_free(H);
        d_free(energies);
        three_vector_incr(kval,kdiff);
    }
    
    g_slice_free(Material, material);
    g_slice_free(Model, model);
        
    if (filename)
        fclose(outstream);

    return 0;
}
