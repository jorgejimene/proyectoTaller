//
// Created by jojimene on 09/04/2025.
//
/**
 * @file  R_z.h
 * @author Jorge
 * @brief Función R_z usada en el proyecto, genera matrices de rotacion al rededor del eje z
 * @date 13-4-2025
 **********************/
#ifndef C___R_Z_H
#define C___R_Z_H
#include "Matrix.h"

/**
 * @details Función que genera una matriz de rotación 3x3 alrededor del eje z, para
 * un ángulo dado en radianes
 * @param angle ángulo de rotación en radianes
 * @return matriz de rotación 3x3
 */
Matrix R_z(double angle);
#endif //C___R_Z_H
