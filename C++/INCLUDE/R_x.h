//
// Created by jojimene on 09/04/2025.
//
/**
 * @file  R_x.h
 * @author Jorge
 * @brief Función R_x usada en el proyecto, genera matrices de rotacion al rededor del eje x
 * @date 13-4-2025
 **********************/

#ifndef C___ROOTMAT_H
#define C___ROOTMAT_H

#include "Matrix.h"

/**
 * @details Función que genera una matriz de rotación 3x3 alrededor del eje x, para
 * un ángulo dado en radianes
 * @param angle angulo de rotación en radianes
 * @return matriz de rotación 3x3
 */
Matrix R_x(double angle);

#endif //C___ROOTMAT_H
