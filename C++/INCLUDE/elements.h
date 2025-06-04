/**
* @file elements.h
 * @brief Cálculo de los elementos keplerianos osculantes a partir del vector de estado.
 */

#ifndef ELEMENTS_H
#define ELEMENTS_H
#include "vector.h"
#include "SAT_Const.h"
#include <math.h>

/**
 * @brief Calcula los elementos orbitales keplerianos osculantes para órbitas elípticas.
 *
 * @param y       Vector de estado del satélite (posición y velocidad) [x, y, z, vx, vy, vz].
 * @param p       [out] Semilatus rectum (m).
 * @param a       [out] Semieje mayor (m).
 * @param e       [out] Excentricidad.
 * @param i       [out] Inclinación orbital (rad).
 * @param Omega   [out] Longitud del nodo ascendente (rad).
 * @param omega   [out] Argumento del periastro (rad).
 * @param M       [out] Anomalía media (rad).
 *
 */
void elements(double y[], double &p, double &a, double &e, double &i,
              double &Omega, double &omega, double &M);

#endif // ELEMENTS_H
