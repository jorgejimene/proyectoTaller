/*
%--------------------------------------------------------------------------
%
% AccelPointMass: Computes the perturbational acceleration due to a point
%				  mass
%
% Inputs:
%   r           Satellite position vector
%   s           Point mass position vector
%   GM          Gravitational coefficient of point mass
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
*/
#include "../INCLUDE/AccelPointMass.h"
#include "cmath"
double norm(const double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void AccelPointMass(double *r,double *s,double GM, double a[3]) {

    double d[3];
    // Relative position vector of satellite w.r.t.point mass
    for (int i = 0; i < 3; i++) {
        d[i] = r[i] - s[i];
    }

    //Acceleration
    for (int i = 0; i < 3; i++) {
        a[i] = -GM * (d[i]/ pow(norm(d),3) + s[i]/pow(norm(s),3));
    }
}