#include <complex.h>
#include "calc-gamma.h"

#ifndef _INC_CALC_THETA_H
#define _INC_CALC_THETA_H

typedef GammaHat ThetaOnePhoton; // same structure

/*typedef struct {
        vector<complex<double> > xx;
        vector<complex<double> > xy;
        vector<complex<double> > xz;
        vector<complex<double> > yx;
        vector<complex<double> > yy;
        vector<complex<double> > yz;
        vector<complex<double> > zx;
        vector<complex<double> > zy;
        vector<complex<double> > zz;
        int Nkc;
        
        ThetaTwoPhoton ( GammaHat ** gh, double omegastep, size_t c, size_t v, int N1, int N2);
        ThetaTwoPhoton (int Nkc); // create one with all zeros
} ThetaTwoPhoton;  // assuming omega_d = 0*/

/*class ThetaTwoPhoton {
    public:
        vector<complex<double> > xx;
        vector<complex<double> > xy;
        vector<complex<double> > xz;
        vector<complex<double> > yx;
        vector<complex<double> > yy;
        vector<complex<double> > yz;
        vector<complex<double> > zx;
        vector<complex<double> > zy;
        vector<complex<double> > zz;
        int Nkc;
        
        ThetaTwoPhoton ( GammaHat ** gh, double omegastep, size_t c, size_t v, int N1, int N2);
        ThetaTwoPhoton (int Nkc); // create one with all zeros
};  // assuming omega_d = 0*/

double calc_theta_get_elapsed_time();

#endif
