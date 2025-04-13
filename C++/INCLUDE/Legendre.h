//
// Created by jojimene on 10/04/2025.
//
/**
 * @file  Legendre.h
 * @author Jorge
 * @brief Función Legendre utilizada en el proyecto, genera una tabla completa de
 * polinomios asociados a Legendre y sus derivadas
 * @date 9-4-2025
 **********************/
#ifndef C___LEGENDRE_H
#define C___LEGENDRE_H
#include <cmath>
#include "Matrix.h"

/**
 * @details Función que calcula los polinomios asociados de Legendre y sus derivadas
 * con respecto al ángulo fi, hasta el grado n y el orden m, guarda los resultados en pnm y dpnm
 * @param n grado máximo de los polinomios a calcular
 * @param m orden máximo
 * @param fi ángulo (en radianes)
 * @param pnm matriz donde se almacenan los polinomios
 * @param dpnm matriz donde se almacenan las derivadas
 */
void Legendre(int n, int m, double fi, Matrix &pnm, Matrix &dpnm);
#endif //C___LEGENDRE_H
