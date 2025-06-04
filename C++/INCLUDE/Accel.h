//
// Created by jojimene on 2025-XX-XX
//

#ifndef C___ACCEL_H
#define C___ACCEL_H
#include "Matrix.h"
#include <SAT_Const.h>
#include "AccelHarmonic.h"
#include "AccelPointMass.h"
#include "AuxParam.h"
#include "GHAMatrix.h"
#include "IERS.h"
#include "JPL_Eph_DE430.h"
#include "Mjday_TBD.h"
#include "NutMatrix.h"
#include "PoleMatrix.h"
#include "PrecMatrix.h"
#include "timediff.h"
/**
 * @brief Calcula la aceleración de un satélite en órbita terrestre debido a:
 *        - el campo de gravedad armónico de la Tierra,
 *        - perturbaciones gravitacionales del Sol y la Luna,
 *        - presión de radiación solar,
 *        - y arrastre atmosférico (si está implementado).
 *
 * @param x Tiempo desde la época (en segundos)
 * @param Y Vector estado del satélite (posición y velocidad) en sistema ICRF/EME2000
 * @return Derivada del vector estado (velocidad y aceleración)
 */

Matrix Accel(double x, const Matrix& Y);

#endif //C___ACCEL_H
