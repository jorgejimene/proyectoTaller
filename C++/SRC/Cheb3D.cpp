#include "../INCLUDE/Cheb3D.h"

#include <iostream>
/*%--------------------------------------------------------------------------
%
% Chebyshev approximation of 3-dimensional vectors
%
% Inputs:
%     N       Number of coefficients
%     Ta      Begin interval
%     Tb      End interval
%     Cx      Coefficients of Chebyshev polyomial (x-coordinate)
%     Cy      Coefficients of Chebyshev polyomial (y-coordinate)
%     Cz      Coefficients of Chebyshev polyomial (z-coordinate)
%
% Last modified:   2018/01/27   M. Mahooti
%
%--------------------------------------------------------------------------
*/
#include <stdexcept>

void Cheb3D(double t, int N, double Ta, double Tb,
            const double Cx[], const double Cy[], const double Cz[],
            double ChebApp[3]) {

    if (t < Ta || t > Tb) {
        throw std::out_of_range("ERROR: Time out of range in Cheb3D::Value");
    }

    double tau = (2.0 * t - Ta - Tb) / (Tb - Ta);

    double f1[3] = {0.0, 0.0, 0.0};
    double f2[3] = {0.0, 0.0, 0.0};
    double old_f1[3];

    for (int i = N - 1; i >= 1; --i) {
        for (int j = 0; j < 3; ++j) {
            old_f1[j] = f1[j];
        }

        f1[0] = 2.0 * tau * f1[0] - f2[0] + Cx[i];
        f1[1] = 2.0 * tau * f1[1] - f2[1] + Cy[i];
        f1[2] = 2.0 * tau * f1[2] - f2[2] + Cz[i];

        for (int j = 0; j < 3; ++j) {
            f2[j] = old_f1[j];
        }
    }

    ChebApp[0] = tau * f1[0] - f2[0] + Cx[0];
    ChebApp[1] = tau * f1[1] - f2[1] + Cy[0];
    ChebApp[2] = tau * f1[2] - f2[2] + Cz[0];
}
