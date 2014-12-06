#include <glib.h>
#include "box.h"
#include "absorption.h"

OnePhoton ** solve_box_simpler (GHashTable * hash, const Box * b, double tol, int max_depth, double max_eps, int num, gboolean bf);
void simple_init();

