/**
* @file EqnEquinox.h
 * @brief Cálculo de la ecuación de los equinoccios.
 */

#ifndef EQNEQUINOX_H
#define EQNEQUINOX_H
#include <valarray>
#include "NutAngles.h"
#include "MeanObliquity.h"
/**
 * @brief Calcula la ecuación de los equinoccios.
 *
 * @param Mjd_TT Fecha juliana modificada (Tiempo Terrestre).
 * @return double Ecuación de los equinoccios (dpsi * cos(eps)) en radianes.
 *
 */
double EqnEquinox(double Mjd_TT);

#endif // EQNEQUINOX_H
