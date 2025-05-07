//
// Created by jojimene on 07/05/2025.
//
#include "Matrix.h"
#include "R_y.h"
#include "R_x.h"

/*
%--------------------------------------------------------------------------
%
% PoleMatrix.m
%
% Purpose:
%   Transformation from pseudo Earth-fixed to Earth-fixed coordinates
%   for a given date
%
% Input:
%   Pole coordinte(xp,yp)
%
% Output:
%   PoleMat   Pole matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
*/
Matrix PoleMatrix (double xp,double yp){
    return R_y(-xp) * R_x(-yp);
}
