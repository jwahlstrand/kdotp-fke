#include "config.h"
#include "three-vector.h"
#include "calc-gamma.h"
#include <gsl/gsl_spline.h>
#include <glib.h>

#ifndef ZENER_H
#define ZENER_H

Gamma ** zener_renormalize (GList * list, double epsilon, ThreeVector * dir);
complex double * zener_renormalize_bundle (GList * list, double epsilon, ThreeVector * dir, Bundle * bundle);
double zener_get_elapsed_time();

#endif
