//
// Created by jojimene on 23/04/2025.
//

#include "../INCLUDE/Unit.h"
#include <cmath>
/*--------------------------------------------------------------------------
%
%  unit.m
%
%  this function calculates a unit vector given the original vector. if a
%  zero vector is input, the vector is set to zero.
%
%  input:
%    vec         - vector
%
%  output:
%    outvec      - unit vector
%
%--------------------------------------------------------------------------
*/
void unit (double vec[], double *outvec) {

    double small = 0.000001;
    double magv = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

    if (magv > small) {
        for (int i = 0; i <= 2; i++) {
            outvec[i] = vec[i] / magv;
        }
    }
    else {
        for (int i = 0; i < 3; i++) {
            outvec[i] = 0.0;
        }
    }
}