#ifndef IERS_H
#define IERS_H

#include "Matrix.h"
#include "Matrix.h"
#include "SAT_Const.h"
/**
 * @brief Interpola parámetros de rotación terrestre a partir de los datos del IERS.
 *
 *
 * @param eop     Matriz con los datos EOP, donde las columnas representan diferentes fechas.
 * @param Mjd_UTC Fecha juliana modificada (UTC).
 * @param x_pole  Coordenada x del polo [rad] (salida).
 * @param y_pole  Coordenada y del polo [rad] (salida).
 * @param UT1_UTC Diferencia UT1 - UTC [s] (salida).
 * @param LOD     Longitud del día [s] (salida).
 * @param dpsi    Nutación en longitud [rad] (salida).
 * @param deps    Nutación en oblicuidad [rad] (salida).
 * @param dx_pole Corrección en x del polo [rad] (salida).
 * @param dy_pole Corrección en y del polo [rad] (salida).
 * @param TAI_UTC Diferencia TAI - UTC [s] (salida).
 * @param interp  Método de interpolación:
 *                - 'l': lineal
 *                - 'n' u otro: sin interpolación (valor directo del día)
 */
void IERS(Matrix &eop,
          double Mjd_UTC,
          double &x_pole,
          double &y_pole,
          double &UT1_UTC,
          double &LOD,
          double &dpsi,
          double &deps,
          double &dx_pole,
          double &dy_pole,
          double &TAI_UTC,
          char interp);

#endif // IERS_H
