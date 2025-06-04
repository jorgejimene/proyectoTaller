//
// Created by jojimene on 03/04/2025.
//
/**
 * @file  NutAngles.h
 * @brief Cálculo de los ángulos de nutación en longitud y oblicuidad
 **********************/

#ifndef C___NUTANGLES_H
#define C___NUTANGLES_H
#include <valarray>
#include "SAT_Const.h"
#include "Matrix.h"
#include <cmath>
/**
 * @details Calcula la nutación en longitud (dpsi) y oblicuidad (deps) para una fecha dada en Tiempo Terrestre (MJD_TT).
 *          Basado en modelos astronómicos de nutación.
 * @param Mjd_TT Fecha en Tiempo Terrestre, expresada en Modified Julian Date (MJD)
 * @param dpsi Referencia a variable donde se almacena la nutación en longitud (radianes)
 * @param deps Referencia a variable donde se almacena la nutación en oblicuidad (radianes)
 */

void NutAngles(double Mjd_TT, double &dpsi, double &deps);

#endif //C___NUTANGLES_H
