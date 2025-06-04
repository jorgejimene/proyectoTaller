/**
* @file EccAnom.h
 * @brief Cálculo de la anomalía excéntrica para órbitas elípticas.
 */

#ifndef ECCANOM_H
#define ECCANOM_H
#include <stdexcept>
#include "cmath"
/**
 * @brief Calcula la anomalía excéntrica a partir de la anomalía media y la excentricidad.
 *
 * @param M  Anomalía media en radianes.
 * @param e  Excentricidad orbital (0 <= e < 1).
 * @return   Anomalía excéntrica en radianes.
 *
 * @throws std::runtime_error si no converge en el número máximo de iteraciones.
 */
double EccAnom(double M, double e);

#endif // ECCANOM_H
