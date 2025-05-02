#include "../INCLUDE/Cheb3D.h"
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
void Cheb3D(double t,int N,double Ta,double Tb,double Cx[],double Cy[],double Cz[], double ChebApp[3]) {
    // Check validity
    if ( (t>Ta) && (Tb>t) ) {
        // Clenshaw algorithm
        double tau = (2*t-Ta-Tb)/(Tb-Ta);

        double f1[] = {0,0,0};
        double f2[] = {0,0,0};

        for (int i=N; i>=1;i--) {
            double old_f1[3];
            for (int j=0; j<3; j++) {
                old_f1[j] = f1[j];
                f1[j] = 2*tau*f1[j]-f2[j];
            }
            f1[0] += Cx[i];
            f1[1] += Cy[i];
            f1[2] += Cz[i];

            for (int j=0; j<3; j++) {
                f2[j] = old_f1[j];
            }
        }
        for (int j=0;j<3;j++) {
            ChebApp[j] = tau*f1[j]-f2[j];
        }
        ChebApp[0] += Cx[0];
        ChebApp[1] += Cy[0];
        ChebApp[2] += Cz[0];
    }
}
