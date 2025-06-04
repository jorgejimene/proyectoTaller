#ifndef G_ACCELHARMONIC_H
#define G_ACCELHARMONIC_H

#include "Matrix.h"  // Asegúrate de que esta cabecera defina correctamente la clase Matrix
#include <iostream>

#include "../INCLUDE/AccelHarmonic.h"
/**
 * @brief Calcula el gradiente del campo gravitacional armónico de la Tierra.
 *
 *
 * @param r     Vector de posición del satélite (en el sistema true-of-date).
 * @param U     Matriz de transformación al sistema fijo del cuerpo (body-fixed).
 * @param n_max Grado máximo del modelo gravitacional.
 * @param m_max Orden máximo del modelo gravitacional.
 * @return      Matriz 3x3 del gradiente de la aceleración respecto a la posición (dA/dr),
 *              en el sistema true-of-date.
 */
Matrix G_AccelHarmonic(Matrix &r, Matrix &U, int n_max, int m_max);

#endif // G_ACCELHARMONIC_H
