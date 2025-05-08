//
// Created by jojimene on 08/05/2025.
//
#include "SAT_Const.h"

/*
%--------------------------------------------------------------------------
%
% Purpose:
%   Computes the osculating Keplerian elements from the satellite state
%   vector for elliptic orbits
%
% Input:
%    y        State vector (x,y,z,vx,vy,vz)
%
% Outputs:
%    p        semilatus rectum [m]
%    a        Semimajor axis
%    e        Eccentricity
%    i        Inclination [rad]
%    Omega    Longitude of the ascending node [rad]
%    omega    Argument of pericenter [rad]
%    M        Mean anomaly [rad]
%
% Notes:
%   The function cannot be used with state vectors describing a circular
%   or non-inclined orbit.
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
#import "cmath"
#include "vector.h"

void elements (double y[],double &p, double &a, double &e, double &i, double &Omega, double &omega,double &M){

        double pi2 = 2*M_PI;
        double r[3]; // Position
        for(int i=0;i<3;i++){
            r[i]=y[i];
        }
        double v[3];                                      // Velocity
        for(int i=3;i<6;i++){
            v[i-3]=y[i];
        }

        double h[3];
        cross(r,v,h);                                    // Areal velocity
        double magh = norm(h);
        p = magh*magh/consts.GM_Earth;
        double H = norm(h);

        Omega = atan2 ( h[0], -h[1] );                     // Long. ascend. node
        Omega = fmod(Omega,pi2);
        i= atan2(sqrt(h[0]*h[0]+h[1]*h[1]), h[2] ); // Inclination
        double u= atan2 ( r[2]*H, -r[0]*h[1]+r[1]*h[1] );    // Arg. of latitude

        double R  = norm(r);                                      // Distance

        a = 1/(2/R-dot(v,v)/consts.GM_Earth);               // Semi-major axis

        double eCosE = 1-R/a;                                     // e*cos(E)
        double eSinE = dot(r,v)/sqrt(consts.GM_Earth*a);           // e*sin(E)

        double e2 = eCosE*eCosE +eSinE*eSinE;
        e  = sqrt(e2);                                     // Eccentricity
        double E  = atan2(eSinE,eCosE);                           // Eccentric anomaly

        M  = fmod(E-eSinE,pi2);                             //Mean anomaly

        double nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          // True anomaly

        omega = fmod(u-nu,pi2);                             // Arg. of perihelion
}
