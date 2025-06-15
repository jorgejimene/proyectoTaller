//
// Created by Jorge on 30/05/2025.
//

#ifndef DEINTEG_H
#define DEINTEG_H
#include <functional>
#include <Matrix.h>
#include <array>
#include <iostream>
#include <Sign_.h>
#include <stdexcept>
#include <limits>
Matrix DEInteg(Matrix (*func)(double x, const Matrix &Y), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y);
#endif //DEINTEG_H
