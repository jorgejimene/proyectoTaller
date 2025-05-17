#include <Matrix.h>
#include <SAT_Const.h>

#include "AuxParam.h"
#include "GHAMatrix.h"
#include "IERS.h"
#include "Mjday_TBD.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "PrecMatrix.h"
#include "timediff.h"
/*
%--------------------------------------------------------------------------
%
% Accel.m
%
% Purpose:
%   Computes the acceleration of an Earth orbiting satellite due to
%    - the Earth's harmonic gravity field,
%    - the gravitational perturbations of the Sun and Moon
%    - the solar radiation pressure and
%    - the atmospheric drag
%
% Inputs:
%   Mjd_TT      Terrestrial Time (Modified Julian Date)
%   Y           Satellite state vector in the ICRF/EME2000 system
%
% Output:
%   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------

void Accel(double x, double Y[6], double Z[6]) {


    Matrix eopdata = Matrix::LoadFromFile("eopdata.txt");
    double x_pole, y_pole,UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    IERS(eopdata, auxParam.Mjd_UTC + x/86400,x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole,TAI_UTC,'l');
    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC, GPS_UTC;
    timediff(UT1_UTC,TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_UT1 = auxParam.Mjd_UTC + x/86400 + UT1_UTC/86400;
    double Mjd_TT = auxParam.Mjd_UTC + x/86400 + TT_UTC/86400;

    Matrix P = PrecMatrix(consts.MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
     r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE430(MJD_TDB);

    % Acceleration due to harmonic gravity field
    a = AccelHarmonic(Y(1:3), E, AuxParam.n, AuxParam.m);

    % Luni-solar perturbations
    if (AuxParam.sun)
        a = a + AccelPointMass(Y(1:3),r_Sun,const.GM_Sun);
    end

    if (AuxParam.moon)
        a = a + AccelPointMass(Y(1:3),r_Moon,const.GM_Moon);
    end

    % Planetary perturbations
    if (AuxParam.planets)
        a = a + AccelPointMass(Y(1:3),r_Mercury,const.GM_Mercury);
    a = a + AccelPointMass(Y(1:3),r_Venus,const.GM_Venus);
    a = a + AccelPointMass(Y(1:3),r_Mars,const.GM_Mars);
    a = a + AccelPointMass(Y(1:3),r_Jupiter,const.GM_Jupiter);
    a = a + AccelPointMass(Y(1:3),r_Saturn,const.GM_Saturn);
    a = a + AccelPointMass(Y(1:3),r_Uranus,const.GM_Uranus);
    a = a + AccelPointMass(Y(1:3),r_Neptune,const.GM_Neptune);
    a = a + AccelPointMass(Y(1:3),r_Pluto,const.GM_Pluto);
    end

    dY = [Y(4:6);a];
}
*/