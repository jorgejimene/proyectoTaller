//
// Created by jojimene on 30/04/2025.
//

#ifndef C___IERS_H
#define C___IERS_H

#include "Matrix.h"

void IERS(Matrix &eop, double Mjd_UTC,  double &x_pole, double &y_pole, double &UT1_UTC, double &LOD, double &dpsi, double &deps, double &dx_pole, double &dy_pole, double &TAI_UTC,char interp = 'o');

#endif //C___IERS_H
