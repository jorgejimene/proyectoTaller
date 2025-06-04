//
// Created by jojimene on 2025-XX-XX
//

#ifndef C___JPL_EPH_DE430_H
#define C___JPL_EPH_DE430_H

#include "Matrix.h"

#include <iostream>

#include "../INCLUDE/Cheb3D.h"
/**
 * @brief Calcula las posiciones ecuatoriales del Sol, la Luna y los nueve planetas
 *        mayores utilizando las efemérides JPL DE430.
 *
 * @param Mjd_TDB Fecha juliana modificada en Tiempo Dinámico Barycéntrico (TDB)
 * @param r_Earth Posición del centro de masa del sistema solar (SSB) - vector 3x1 [m]
 * @param r_Mars Posición de Marte - vector 3x1 [m]
 * @param r_Mercury Posición de Mercurio - vector 3x1 [m]
 * @param r_Venus Posición de Venus - vector 3x1 [m]
 * @param r_Jupiter Posición de Júpiter - vector 3x1 [m]
 * @param r_Saturn Posición de Saturno - vector 3x1 [m]
 * @param r_Uranus Posición de Urano - vector 3x1 [m]
 * @param r_Neptune Posición de Neptuno - vector 3x1 [m]
 * @param r_Pluto Posición de Plutón - vector 3x1 [m]
 * @param r_Moon Posición de la Luna - vector 3x1 [m]
 * @param r_Sun Posición del Sol (geocéntrica, ICRF) - vector 3x1 [m]
 */
void JPL_Eph_DE430(double Mjd_TDB, Matrix& r_Earth, Matrix& r_Mars, Matrix& r_Mercury, Matrix& r_Venus,
                   Matrix& r_Jupiter, Matrix& r_Saturn, Matrix& r_Uranus, Matrix& r_Neptune,
                   Matrix& r_Pluto, Matrix& r_Moon, Matrix& r_Sun);

#endif //C___JPL_EPH_DE430_H
