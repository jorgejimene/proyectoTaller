#ifndef GHAMATRIX_H
#define GHAMATRIX_H

#include "../INCLUDE/matrix.h"  // Asegúrate de incluir la definición de Matrix
#include <../INCLUDE/R_z.h>

#include "../INCLUDE/gast.h"
/**
 * @brief Calcula la matriz de ángulo horario de Greenwich (GHA).
 *
 *
 * @param Mjd_UT1 Fecha Juliana Modificada en UT1.
 * @return Matriz de ángulo horario de Greenwich (3x3).
 */
Matrix GHAMatrix(double Mjd_UT1);

#endif // GHAMATRIX_H
