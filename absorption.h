#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include "tensors.h"
#include "three-vector.h"
#include "calc-gamma.h"
#include "calc-theta.h"
#include <glib.h>

#ifndef ABSORPTION_H
#define ABSORPTION_H

typedef struct {
  int N1, N2;
  double e; // epsilon
  double *o;
  Tensor2 * data;
} OnePhoton;

void one_photon_scale(OnePhoton *o, double weight);
void one_photon_add_and_scale(OnePhoton *o1, OnePhoton *o2, double weight);
void one_photon_write(OnePhoton *o, FILE *file);
void one_photon_free (OnePhoton *op);

OnePhoton ** blank_abs (const double max_eps, const int num, gboolean include_zero);
OnePhoton ** calculate_absorption (GList * list, const double max_eps, const int num, ThreeVector * dir, gboolean include_zero);

#endif
