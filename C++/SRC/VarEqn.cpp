#include "../INCLUDE/VarEqn.h"

#include <../INCLUDE/SAT_Const.h>

#include "GHAMatrix.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "../INCLUDE/PrecMatrix.h"
#include "../INCLUDE/timediff.h"
#include "../INCLUDE/IERS.h"
#include "../INCLUDE/AuxParam.h"
/*%------------------------------------------------------------------------------
%
% VarEqn.m
%
% Purpose:
%   Computes the variational equations, i.e. the derivative of the state vector
%   and the state transition matrix
%
% Input:
%   x           Time since epoch in [s]
%   yPhi        (6+36)-dim vector comprising the state vector (y) and the
%               state transition matrix (Phi) in column wise storage order
%
% Output:
%   yPhip       Derivative of yPhi
%
% Last modified:   2015/08/12   M. Mahooti
%
%------------------------------------------------------------------------------
*/
Matrix VarEqn(double x, Matrix yPhi) {
    Matrix mat = Matrix::LoadFromFile("eop19620101.txt");
    Matrix matTras = Matrix(mat.getCol(), mat.getFil());
    matTras = mat.transpose();
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, TT_GPS, GPS_TAI,GPS_UTC;

    IERS(matTras, 37668, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, 'l');

    timediff(UT1_UTC,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
    double Mjd_UT1 = auxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    // Transformation matrix
    Matrix P = PrecMatrix(consts.MJD_J2000,auxParam.Mjd_TT + x/86400);
    Matrix N = NutMatrix(auxParam.Mjd_TT + x/86400);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    Matrix r(3,1), v(3,1);
    for (int i=1;i<=3;i++) {
        r(i,1) = yPhi(i,1);
    }
    for (int i=4;i<=6; i++) {
        v(i,1) = yPhi(i,1);
    }
    Matrix Phi(6,6);

    // State transition matrix
    for (int j=1;j<=6;j++) {
        Matrix col(6,1);
        for (int i=1; i<=6; i++) {

        }
    }

    % Acceleration and gradient
    a = AccelHarmonic ( r, E, AuxParam.n, AuxParam.m );
    G = G_AccelHarmonic ( r, E, AuxParam.n, AuxParam.m );

    % Time derivative of state transition matrix
    yPhip = zeros(42,1);
    dfdy = zeros(6);

    for i=1:3
        for j=1:3
            dfdy(i,j) = 0.0;                 % dv/dr(i,j)
            dfdy(i+3,j) = G(i,j);            % da/dr(i,j)
            if ( i==j )
                dfdy(i,j+3) = 1;
            else
                dfdy(i,j+3) = 0;             % dv/dv(i,j)
            end
            dfdy(i+3,j+3) = 0.0;             % da/dv(i,j)
        end
    end

    Phip = dfdy*Phi;

    % Derivative of combined state vector and state transition matrix
    for i=1:3
        yPhip(i)   = v(i);                 % dr/dt(i)
        yPhip(i+3) = a(i);                 % dv/dt(i)
    end

    for i=1:6
        for j=1:6
            yPhip(6*j+i) = Phip(i,j);     % dPhi/dt(i,j)
        end
    end
}
