#include "k-dot-p.h"

#ifndef ZB_H
#define ZB_H

Model * zincblende2pba_new ();

complex double * Px8 ();
complex double * Py8 ();
complex double * Pz8 ();

Model * zincblende8sph_new ();
Model * zincblende8sph2_new ();
Model * zincblende8_new ();

complex double * Px14 ();
complex double * Py14 ();
complex double * Pz14 ();

void px14 (complex double *);

Model * zincblende14ns_new ();
Model * zincblende14nr_new ();
Model * zincblende14_new ();
Model * zincblende30_new ();

void bundle8 (Model *);
void bundle14 (Model *);

#endif
