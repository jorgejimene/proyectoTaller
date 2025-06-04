//
// Created by jojimene on 07/05/2025.
//
/**
 * @file  PoleMatrix.h
 * @brief Calcula la matriz de transformación de coordenadas de pseudo Earth-fixed a Earth-fixed para una fecha dada
 **********************/

#ifndef C___POLEMATRIX_H
#define C___POLEMATRIX_H

#include "Matrix.h"
#include "Matrix.h"
#include "R_y.h"
#include "R_x.h"
/**
 * @details Calcula la matriz de polo que transforma coordenadas pseudo Earth-fixed a Earth-fixed
 *          utilizando los ángulos de polo xp e yp.
 * @param xp Ángulo de polo en la coordenada x (radianes)
 * @param yp Ángulo de polo en la coordenada y (radianes)
 * @return Matriz de transformación PoleMatrix
 */

Matrix PoleMatrix(double xp, double yp);

#endif //C___POLEMATRIX_H
