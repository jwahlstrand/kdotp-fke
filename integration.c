#include <stdlib.h>
#include <glib/gprintf.h>
#include <math.h>
#include "matrix-element-pert.h"
#include "absorption.h"
#include "integration.h"
#include "zener.h"

extern Model * model;

//#define BAD_BOX_CHECK

int depth;
double fraction_done;
Box * integration_box;
gboolean include_zero;

/* following two functions associate a kperp with a string */
/* and vice versa */

static
gchar * kperp2string (const ThreeVector * kp) {
	gchar buffer[40];
	g_sprintf(buffer, "%1.6e_%1.6e", kp->x[0], kp->x[1]);
	return g_strdup(buffer);
}

/* find the absorption spectrum for one kperp */
/* look up kperp in the hash table first */
/* if it's not there, calculate it, then add it to the hash table */

OnePhoton ** solve_kperp (GHashTable * hash, ThreeVector * kp, ThreeVector * normal, double max_eps, int num) {
	gchar* key;
	OnePhoton ** solution;
	
    /*
    #ifdef BAD_BOX_CHECK    // need to redo this, restriction is not as great once we get 
    	g_assert(fabs(fabs(kp->x[0])-fabs(kp->x[1]))>1e-6);  // bad box filtering should prevent this from happening
    	g_assert(fabs(kp->x[0])>1e-6);
    	g_assert(fabs(kp->x[1])>1e-6);
    #endif
    */
    
    /*if(fabs(fabs(kp->x[0])-fabs(kp->x[1]))<1e-6) {
	    solution = blank_abs (max_eps, num, calc_type);
	    return solution;
	}*/

    key = kperp2string (kp);
    gpointer p1 = NULL;
    if (hash) {
        p1 = g_hash_table_lookup (hash, key);
    }
    if (!hash || p1 == NULL) { /* spectrum wasn't found in hashtable, so find it... */
        GList * list;
        g_message("calculating spectrum: %e,%e",kp->x[0],kp->x[1]);
        list = calc_w_phi_coeffs (kp, normal);
        
        if (list==NULL) {  // we have a degeneracy problem
            solution = blank_abs (max_eps, num, include_zero);
        }
        else {
            solution = calculate_absorption (list, max_eps, num, normal, include_zero);
            matrix_element_list_free (list);
        }
        
        if (hash) 
            g_hash_table_insert (hash, key, solution);  /* ...and add to hash */
    }
    else {            /* spectrum was found, so use it */
        g_message("found spectrum: %e,%e",kp->x[0],kp->x[1]);
        solution = (OnePhoton **) p1;
        g_free(key);  /* since we won't be putting it in the hash */
    }
    
    return solution;
}

/* GHashTable funcs */

static
void destroy_value (gpointer value) {
	g_assert(value);
	OnePhoton * a = (OnePhoton *) value;
	//one_photon_free(a);
}

/* * * * * * * * * * * * * * * * * * */

/* used to clean the hash table after solving a box */

void simple_init() {
    depth=0;
    fraction_done=0.0;
}

OnePhoton ** solve_box_simpler (GHashTable * hash, const Box * b, double tol, int max_depth, double max_eps, int num, gboolean incl_zero) {
    OnePhoton **asum = NULL;
    OnePhoton **a1;
    
    include_zero = incl_zero;

    /* divide box into 4**max_depth boxes, sum over them */
    
    double kx, ky;
    
    const double kxmin = (b->center->x[0])-(b->width/2);
    const double kxmax = (b->center->x[0])+(b->width/2);
    const double kymin = (b->center->x[1])-(b->width/2);
    const double kymax = (b->center->x[1])+(b->width/2);
    
    if(max_depth==0)
      return solve_kperp (NULL, b->center, b->normal, max_eps, num);
    
    const double width = (b->width)/pow(2.0,max_depth);
    int d;
    for (kx = kxmin+width/2;kx<kxmax;kx=kx+width) {
        for (ky = kymin+width/2;ky<kymax;ky=ky+width) {
            ThreeVector * c = three_vector_new (kx, ky, 0);
            Box * box1 = new_box (c, b->normal, width);
            if (!asum) {
                asum=solve_kperp (NULL, box1->center, box1->normal, max_eps, num);
                for(d=0;d<num;d++) {
                  one_photon_scale(asum[d],width*width);
                }
            }
            else {
                a1=solve_kperp (NULL, box1->center, box1->normal, max_eps, num);
                for(d=0;d<num;d++)
                  one_photon_add_and_scale(asum[d],a1[d], width*width);
                destroy_value(a1);
            }
            g_slice_free(Box, box1);
	        
	        fraction_done += pow(1/4.0,max_depth);  // only increment if this is the deepest level
        }
    }

    return asum;
}
