//
// Created by jojimene on 09/04/2025.
//
/**
 * @file  R_y.h
 * @author Jorge
 * @brief Función R_y usada en el proyecto, genera matrices de rotacion al rededor del eje y
 * @date 13-4-2025
 **********************/
#ifndef C___R_Y_H
#define C___R_Y_H

#include "Matrix.h"
#include <cmath>
/**
 * @details Función que genera una matriz de rotación 3x3 alrededor del eje y, para
 * un ángulo dado en radianes
 * @param angle ángulo de rotación en radianes
 * @return matriz de rotación 3x3
 */
Matrix R_y(double angle);
#endif //C___R_Y_H
