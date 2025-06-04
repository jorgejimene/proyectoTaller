//
// Creado por jojimene el 07/05/2025.
//

#ifndef C___NUTMATRIX_H
#define C___NUTMATRIX_H

#include "Matrix.h"
#include "MeanObliquity.h"
#include "NutAngles.h"
#include "R_x.h"
#include "R_z.h"

/**
 * @brief Calcula la matriz de nutación.
 *
 * @param Mjd_TT Fecha Juliana Modificada en Tiempo Terrestre (TT)
 * @return Matriz de transformación de nutación (3x3)
 */
Matrix NutMatrix(double Mjd_TT);

#endif //C___NUTMATRIX_H
