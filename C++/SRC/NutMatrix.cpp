//
// Created by jojimene on 07/05/2025.
//
/*
%--------------------------------------------------------------------------
%
% NutMatrix.m
%
% Purpose:
%   Transformation from mean to true equator and equinox
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%   NutMat    Nutation matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/

#include "NutMatrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"

Matrix NutMatrix (double Mjd_TT){
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity(Mjd_TT);

    // Nutation in longitude and obliquity
    double dpsi,deps;
    NutAngles (Mjd_TT, dpsi, deps);

    // Transformation from mean to true equator and equinox

    Matrix NutMat = R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);

    return NutMat;
}
