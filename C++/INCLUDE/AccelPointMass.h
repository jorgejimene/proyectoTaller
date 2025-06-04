/**
* @file AccelPointMass.h
 * @brief Cálculo de la aceleración perturbacional debida a una masa puntual.
 */

#ifndef ACCELPOINTMASS_H
#define ACCELPOINTMASS_H

#include "Matrix.h"

/**
 * @brief Calcula la aceleración perturbacional debida a una masa puntual.
 *
 * La aceleración se calcula como la diferencia entre el vector de posición del satélite y el de la masa perturbadora.
 *
 * @param r Vector de posición del satélite (3x1).
 * @param s Vector de posición de la masa puntual perturbadora (3x1).
 * @param GM Coeficiente gravitacional de la masa perturbadora [m³/s²].
 * @return Aceleración perturbacional (3x1).
 */
Matrix AccelPointMass(Matrix &r, Matrix &s, double GM);

#endif // ACCELPOINTMASS_H
