#include <complex.h>
#include <gsl/gsl_odeiv.h>

#ifndef PERTODE_H
#define PERTODE_H

typedef struct {
    const gsl_odeiv_step_type * Te;
    gsl_odeiv_step * se;
    gsl_odeiv_control * ce;
    gsl_odeiv_evolve * ee;
    gsl_odeiv_system syse;
} pert_ode;

pert_ode * pert_ode_init (size_t Nv);
void pert_ode_reset (pert_ode * po);
int pert_ode_evolve (pert_ode * po, double * te, double dk, double * he, double data[]);
void pert_ode_free (pert_ode * po);

double * pert_initial_condition (complex double * mat, size_t Nv);

#endif
