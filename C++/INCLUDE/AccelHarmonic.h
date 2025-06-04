/**
* @file AccelHarmonic.h
 * @brief Cálculo de la aceleración debida al campo gravitatorio armónico de un cuerpo central.
 */

#ifndef ACCELHARMONIC_H
#define ACCELHARMONIC_H

#include "Matrix.h"
#include <Legendre.h>

#include "globalMatrix.h"

/**
 * @brief Calcula la aceleración gravitatoria debida al campo armónico de un cuerpo central.
 *
 * Esta función utiliza coeficientes armónicos (Cnm, Snm) hasta un grado y orden máximos
 * para calcular la aceleración gravitatoria en un sistema inercial.
 *
 * @param r Vector de posición del satélite en el sistema inercial (3x1).
 * @param E Matriz de transformación del sistema inercial al sistema fijo del cuerpo (3x3).
 * @param n_max Grado máximo del desarrollo armónico.
 * @param m_max Orden máximo del desarrollo armónico.
 * @return Aceleración del satélite (3x1) en el sistema inercial.
 */
Matrix AccelHarmonic(Matrix &r, Matrix &E, int n_max, int m_max);

#endif // ACCELHARMONIC_H
