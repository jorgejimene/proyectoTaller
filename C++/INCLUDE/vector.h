//
// Created by jojimene on 08/05/2025.
//
/**
 * @file  vector.h
 * @brief Funciones para operaciones vectoriales: norma, producto punto y producto cruz
 **********************/

#ifndef C___VECTOR_H
#define C___VECTOR_H
#include <cmath>
/**
 * @brief Calcula la norma (módulo) de un vector de 3 componentes.
 * @param v Vector de entrada (3 elementos)
 * @return Norma del vector
 */
double norm(const double v[3]);

/**
 * @brief Calcula el producto punto de dos vectores de 3 componentes.
 * @param r Primer vector (3 elementos)
 * @param x Segundo vector (3 elementos)
 * @return Producto punto (escalar)
 */
double dot(const double r[], const double x[]);

/**
 * @brief Calcula el producto cruz de dos vectores de 3 componentes.
 * @param v1 Primer vector (3 elementos)
 * @param v2 Segundo vector (3 elementos)
 * @param v Vector resultado (3 elementos)
 */
void cross(const double v1[], const double v2[], double v[]);

#endif //C___VECTOR_H
