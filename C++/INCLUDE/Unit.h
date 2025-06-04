//
// Created by jojimene on 23/04/2025.
//
/**
 * @file  Unit.h
 * @brief Función unit que calcula el vector unitario a partir de un vector dado.
 **********************/

#ifndef C___UNIT_H
#define C___UNIT_H
#include <cmath>
/**
 * @brief Calcula el vector unitario de un vector dado. Si el vector de entrada es cero, el resultado es un vector cero.
 *
 * @param vec Vector de entrada (arreglo de 3 elementos)
 * @param outvec Vector de salida (arreglo de 3 elementos) donde se almacena el vector unitario
 */
void unit(double vec[], double *outvec);

#endif //C___UNIT_H
