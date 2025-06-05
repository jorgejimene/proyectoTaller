
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
%---------------------------------------------------------------------------
*/

#include <../INCLUDE/Accel.h>
using namespace std;
Matrix Accel(double x, const Matrix& Y) {
    Matrix eopdata = Matrix::LoadFromFile("eop19620101.txt");
    Matrix eopTrans = eopdata.transpose();
    double x_pole, y_pole,UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    IERS(eopTrans, auxParam.Mjd_UTC + x/86400,x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole,TAI_UTC,'l');
    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC, GPS_UTC;
    cout << "accel" << endl;

    timediff(UT1_UTC,TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_UT1 = auxParam.Mjd_UTC + x/86400 + UT1_UTC/86400;
    double Mjd_TT = auxParam.Mjd_UTC + x/86400 + TT_UTC/86400;
    Matrix P = PrecMatrix(consts.MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
    cout << "accel" << endl;

    Matrix r_Earth(3,1);
    Matrix r_Mars(3,1);
    Matrix r_Mercury(3,1);
    Matrix r_Venus(3,1);
    Matrix r_Jupiter(3,1);
    Matrix r_Saturn(3,1);
    Matrix r_Uranus(3,1);
    Matrix r_Neptune(3,1);
    Matrix r_Pluto(3,1);
    Matrix r_Moon(3,1);
    Matrix r_Sun(3,1);

    JPL_Eph_DE430(MJD_TDB, r_Earth,r_Mars,r_Mercury,r_Venus,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun);
    // Acceleration due to harmonic gravity field
    cout << "accel" << endl;
    cout << "dimensiones Y" << endl;
    cout << Y.getFil() << Y.getCol() << endl;
    Matrix r(3,1);
    for (int i = 1; i <= 3; ++i) {
        r(i,1) = Y(i,1);
    }
    r.print();
    cout << "accel" << endl;

    Matrix a = AccelHarmonic(r, E, auxParam.n, auxParam.m);
    // Luni-solar perturbations
    if (auxParam.sun) {
        a = a + AccelPointMass(r, r_Sun,consts.GM_Sun);
    }
    if (auxParam.moon) {
        a = a + AccelPointMass(r,r_Moon,consts.GM_Moon);
    }

    cout << "accel" << endl;
    // Planetary perturbations
    if (auxParam.planets){
        a = a + AccelPointMass(r,r_Mercury,consts.GM_Mercury);
        a = a + AccelPointMass(r,r_Venus,consts.GM_Venus);
        a = a + AccelPointMass(r,r_Mars,consts.GM_Mars);
        a = a + AccelPointMass(r,r_Jupiter,consts.GM_Jupiter);
        a = a + AccelPointMass(r,r_Saturn,consts.GM_Saturn);
        a = a + AccelPointMass(r,r_Uranus,consts.GM_Uranus);
        a = a + AccelPointMass(r,r_Neptune,consts.GM_Neptune);
        a = a + AccelPointMass(r,r_Pluto,consts.GM_Pluto);
    }
    a.print();
    Matrix dY(6,1);
    for (int i=1;i<=3;++i) {
        dY(i,1) = Y(i+3,1);
        dY(i+3,1) = a(i,1);
    }
    return dY;
}
