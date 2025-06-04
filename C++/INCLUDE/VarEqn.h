//
// Created by jojimene on 04/06/2025.
//
/**
 * @file  VarEqn.h
 * @brief Cálculo de las ecuaciones variacionales para derivadas del vector estado y la matriz de transición de estado.
 **********************/

#ifndef C___VAREQN_H
#define C___VAREQN_H

#include "Matrix.h"
#include "AccelHarmonic.h"
#include "GHAMatrix.h"
#include "G_AccelHarmonic.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "../INCLUDE/PrecMatrix.h"
#include "../INCLUDE/timediff.h"
#include "../INCLUDE/IERS.h"
#include "../INCLUDE/AuxParam.h"
/**
 * @brief Calcula la derivada del vector de estado y la matriz de transición del estado.
 *
 * @param x Tiempo desde la época en segundos
 * @param yPhi Vector (6+36) que contiene el vector de estado (6) y la matriz de transición de estado (36) en orden columna
 * @return Matriz que contiene la derivada del vector combinado yPhi
 */
Matrix VarEqn(double x, const Matrix &yPhi);

#endif //C___VAREQN_H
