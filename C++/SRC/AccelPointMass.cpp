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

Matrix AccelPointMass(Matrix &r, Matrix &s, double GM) {
    // Relative position vector d = r - s
    Matrix d = r - s;
    double norm_d = d.normVec();
    double norm_s = s.normVec();
    double norm_d_cubed = norm_d * norm_d * norm_d;
    double norm_s_cubed = norm_s * norm_s * norm_s;

    Matrix term1(3,1), term2(3,1);
    for (int i = 1; i <= 3; ++i) {
        term1(i, 1) = d(i, 1) / norm_d_cubed;
        term2(i, 1) = s(i, 1) / norm_s_cubed;
    }

    Matrix a = term1 + term2;
    for (int i = 1; i <= 3; ++i) {
        a(i, 1) *= -GM;
    }

    return a;
}