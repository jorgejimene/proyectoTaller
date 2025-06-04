/*
*%--------------------------------------------------------------------------
%
%  Fractional part of a number (y=x-[x])
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
 */
#include "../INCLUDE/Frac.h"
double Frac(double x){
    int x2 = (int) x;
    return x-x2;
}