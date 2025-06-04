/**
* @file MeasUpdate.h
 * @brief Actualización de estado y covarianza mediante filtro de Kalman.
 *
 * @details
 * Esta función realiza la actualización de la estimación del estado y la covarianza
 * en un filtro de Kalman, a partir de una medición nueva.
 *
 * @param[in,out] x Vector de estado a actualizar.
 * @param[in] z Vector de medición.
 * @param[in] g Vector de predicción de medición.
 * @param[in] s Vector de desviaciones estándar de la medición.
 * @param[in] G Matriz de sensibilidad (Jacobiana) de la medición respecto al estado.
 * @param[in,out] P Matriz de covarianza del estado a actualizar.
 * @param[in] n Dimensión del vector de estado.
 * @param[out] K Matriz de ganancia de Kalman calculada.
 */

#ifndef MEASUPDATE_H
#define MEASUPDATE_H

#include "Matrix.h"  // Asumiendo que la clase Matrix está declarada aquí

void MeasUpdate(Matrix &x, Matrix &z, Matrix &g, Matrix &s,
                Matrix &G, Matrix &P, int n, Matrix &K);

#endif // MEASUPDATE_H
