/*--------------------------------------------------------------------------
%
% G_AccelHarmonic.m
%
% Purpose:
%   Computes the gradient of the Earth's harmonic gravity field
%
% Inputs:
%   r           Satellite position vector in the true-of-date system
%   U           Transformation matrix to body-fixed system
%   n           Gravity model degree
%   m 			Gravity model order
%
% Output:
%   G    		Gradient (G=da/dr) in the true-of-date system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
#include "../INCLUDE/G_AccelHarmonic.h"

#include <iostream>

#include "../INCLUDE/AccelHarmonic.h"

Matrix G_AccelHarmonic(Matrix &r, Matrix &U, int n_max, int m_max) {
    double d = 1.0;   // Position increment [m]

    Matrix G(3,3);
    Matrix dr (3,1);

    // Gradient
    for (int i=1;i<=3;++i) {
        for (int j=1;j<=3;++j) {
            dr(j,1) = 0.0;
        }
        dr(i,1) = d;
        // Acceleration difference
        Matrix r_plus(3,1), r_minus(3,1);
        for (int j=1;j<=3;++j) {
            r_plus(j,1) = r(j,1) + 0.5 *dr(j,1);
            r_minus(j,1) = r(j,1) - 0.5 *dr(j,1);
        }
        Matrix a_plus = AccelHarmonic(r_plus,U,n_max,m_max);
        Matrix a_minus = AccelHarmonic(r_minus,U,n_max,m_max);
        // Derivative with respect to i-th axis

        for (int j=1;j<=3;++j) {
            G(j,i) = (a_plus(j,1) - a_minus(j,1))/d;
        }
    }

    return G;
}