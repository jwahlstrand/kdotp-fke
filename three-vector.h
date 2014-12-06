#ifndef THREEVECT_H
#define THREEVECT_H

/* represents a three-vector with double precision */
/* this is a little silly */

typedef struct {
  double x[3];
} ThreeVector;

/* represents a three-vector with double precision */

void three_vector_free(ThreeVector *tv);

ThreeVector * three_vector_new (double x, double y, double z);
ThreeVector * three_vector_copy (const ThreeVector *);
ThreeVector * three_vector_add (const ThreeVector * k1, const ThreeVector * k2);
ThreeVector * three_vector_diff (const ThreeVector * k1, const ThreeVector * k2);
void three_vector_incr (ThreeVector *, const ThreeVector *dk);
void three_vector_scale (ThreeVector *, double );
ThreeVector * three_vector_from_string (char * text);
double three_vector_length (const ThreeVector *tv);
ThreeVector * dir_from_string (char * string);
ThreeVector * three_vector_unit (const ThreeVector *tv);
ThreeVector * three_vector_from_string_array (char ** toks);

ThreeVector * three_vector_e (const ThreeVector *g);
ThreeVector * three_vector_f (const ThreeVector *g, const ThreeVector *e );

int trajectory_intersects_bad (const ThreeVector * kperp, const ThreeVector * g);

void three_vector_array_project_inplace (double *P, const double * x, const double * y, const double * z, const ThreeVector * dir, size_t N);
double * three_vector_array_project (const double * x, const double * y, const double * z, const ThreeVector * k, size_t N);

#endif
