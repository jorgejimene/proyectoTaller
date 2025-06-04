/**
* @file LTC.h
 * @brief Transformación del sistema del meridiano de Greenwich a coordenadas tangentes locales.
 *
 * @details Esta función genera una matriz de rotación que transforma coordenadas
 * geodésicas (longitud y latitud) al sistema de coordenadas tangente local
 * (Este-Norte-Zenit, ENZ).
 *
 * Basado en la implementación original en MATLAB por M. Mahooti (2015/08/12).
 */

#ifndef LTC_H
#define LTC_H

#include "Matrix.h"  ///< Definición de la clase Matrix
#include <R_y.h>
#include <R_z.h>
/**
 * @brief Calcula la matriz de transformación a coordenadas locales tangentes.
 *
 * @param lon Longitud geodésica este en radianes.
 * @param lat Latitud geodésica en radianes.
 * @return Matrix Matriz de rotación de 3x3 que transforma del sistema
 * del ecuador terrestre y meridiano de Greenwich al sistema local tangente (Este-Norte-Zenit).
 *
 */
Matrix LTC(double lon, double lat);

#endif // LTC_H
