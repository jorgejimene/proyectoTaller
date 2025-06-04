//
// Created by jojimene on 07/05/2025.
//
/**
 * @file  PrecMatrix.h
 * @brief Matriz de transformación de precesión para coordenadas ecuatoriales
 **********************/

#ifndef C___PRECMATRIX_H
#define C___PRECMATRIX_H
#include "../INCLUDE/SAT_Const.h"
#include "../INCLUDE/R_z.h"
#include "../INCLUDE/R_y.h"
/**
 * @details Calcula la matriz de precesión que transforma coordenadas ecuatoriales
 *          de una época dada (MJD_1) a otra (MJD_2), usando fechas julianas modificadas TT.
 * @param Mjd_1 Época inicial (Modified Julian Date TT)
 * @param Mjd_2 Época final a la que se precesa (Modified Julian Date TT)
 * @return Matriz de transformación de precesión (Matrix)
 */

Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif //C___PRECMATRIX_H
