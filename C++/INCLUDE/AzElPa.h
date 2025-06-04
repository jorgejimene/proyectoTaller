/**
* @file AzElPa.h
 * @brief Cálculo de acimut, elevación y derivadas parciales desde coordenadas topocéntricas.
 */

#ifndef AZELPA_H
#define AZELPA_H
#include "cmath"
/**
 * @brief Calcula el acimut (Az), la elevación (El) y sus derivadas parciales respecto al vector s.
 *
 * @param s Vector de entrada en coordenadas topocéntricas locales (Este-Norte-Zenit) [3 elementos].
 * @param Az Acimut de salida [rad].
 * @param El Elevación de salida [rad].
 * @param dAds Derivadas parciales del acimut respecto a s [3 elementos].
 * @param dEds Derivadas parciales de la elevación respecto a s [3 elementos].
 */
void AzElPa(double s[3], double &Az, double &El, double dAds[3], double dEds[3]);

#endif // AZELPA_H
