//
// Created by jojimene on 07/05/2025.
//
/*
%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the eccentric anomaly for elliptic orbits
%
% Inputs:
%   M         Mean anomaly in [rad]
%   e         Eccentricity of the orbit [0,1]
%
% Output:
%             Eccentric anomaly in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
#include "../INCLUDE/EccAnom.h"

double EccAnom (double M,double e) {

    int maxit = 15;
    int i = 1;
    double E;

    // Starting value
    double pi = 2 * M_PI;
    M = fmod(M, pi);

    if (e < 0.8) {
        E = M;
    } else {
        E = M_PI;
    }
    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    //Iteration
    while (fabs(f) > 1e2 * 2.22044604925031e-16) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i = i + 1;
        if (i == maxit) {
            throw std::runtime_error("Convergence problems in EccAnom");
        }
    }
    return E;
}