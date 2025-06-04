#include "../INCLUDE/AccelHarmonic.h"

#include <cmath>
#include <iostream>

/*--------------------------------------------------------------------------
%
% AccelHarmonic.m
%
% Purpose:
%   Computes the acceleration due to the harmonic gravity field of the
%   central body
%
% Inputs:
%   r           Satellite position vector in the inertial system
%   E           Transformation matrix to body-fixed system
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%
% Output:
%   a           Acceleration (a=d^2r/dt^2)
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
using namespace std;
Matrix AccelHarmonic(Matrix &r, Matrix &E, int n_max, int m_max) {
    double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
    double gm    = 398600.4415e9; // [m^3/s^2]; GGM03S

    // Body-fixed position
    Matrix r_bf = E * r;

    // Auxiliary quantities
    double d = r_bf.norm();// distance
    double latgc = asin(r_bf(3,1)/d);
    double lon = atan2(r_bf(2,1),r_bf(1,1));

    Matrix pnm(n_max+1,m_max+1), dpnm(n_max+1,m_max+1);
    Legendre(n_max,m_max,latgc,pnm,dpnm);


    double dUdr = 0;
    double dUdlatgc = 0;
    double dUdlon = 0;
    double q3 = 0, q2 = q3, q1 = q2;
    double b1, b2, b3;
    for (int n=0;n<=n_max;++n) {
        b1 = (-gm/(d*d))*pow((r_ref/d),n)*(n+1);
        b2 =  (gm/d)*pow(r_ref/d,n);
        b3 =  (gm/d)*pow(r_ref/d,n);
        for (int m=0;m<=std::min(m_max,n);++m) {
            q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
            q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
            q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
        }
        dUdr     = dUdr     + q1*b1;
        dUdlatgc = dUdlatgc + q2*b2;
        dUdlon   = dUdlon   + q3*b3;
        q3 = 0; q2 = q3; q1 = q2;
    }

    // Body-fixed acceleration
    double r2xy = pow(r_bf(1,1),2)+pow(r_bf(2,1),2);

    double ax = (1/d*dUdr-r_bf(3,1)/(d*d*sqrt(r2xy))*dUdlatgc)*r_bf(1,1)-(1/r2xy*dUdlon)*r_bf(2,1);
    double ay = (1/d*dUdr-r_bf(3,1)/(d*d*sqrt(r2xy))*dUdlatgc)*r_bf(2,1)+(1/r2xy*dUdlon)*r_bf(1,1);
    double az =  1/d*dUdr*r_bf(3,1)+sqrt(r2xy)/(d*d)*dUdlatgc;
    Matrix a_bf(3,1);
    a_bf(1,1)=ax;
    a_bf(2,1)=ay;
    a_bf(3,1)=az;
    
    Matrix Etrans = E.transpose();
    // Inertial acceleration
    a_bf= Etrans * a_bf;
    return a_bf;
}
