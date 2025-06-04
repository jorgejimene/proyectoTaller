/**
* @file MeanObliquity.h
 * @brief Cálculo de la oblicuidad media de la eclíptica.
 *
 */

#ifndef MEANOBLIQUITY_H
#define MEANOBLIQUITY_H
#include "SAT_Const.h"

/**
 * @brief Calcula la oblicuidad media de la eclíptica.
 *
 * @param Mjd_TT Fecha Juliana Modificada en Tiempo Terrestre (TT).
 * @return double Oblicuidad media de la eclíptica en radianes.
 *
 * @note Utiliza la constante MJD_J2000 para referencia al 1 de enero de 2000 TT.
 */
double MeanObliquity(double Mjd_TT);

#endif // MEANOBLIQUITY_H
