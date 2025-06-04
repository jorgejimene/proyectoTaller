//
// Created by jojimene on 10/04/2025.
//
/**
 * @file  TimeUpdate.h
 * @brief Función TimeUpdate para actualizar matrices de covarianza y transición
 **********************/

#ifndef C___TIMEUPDATE_H
#define C___TIMEUPDATE_H

#include "Matrix.h"

/**
 * @brief Actualiza la matriz de covarianza P usando la matriz de transición Phi y la matriz Qdt opcional
 *
 * @param P Puntero a la matriz de covarianza que será actualizada
 * @param Phi Puntero a la matriz de transición
 * @param Qdt Puntero a la matriz de ruido Q multiplicado por dt (opcional, puede ser nullptr)
 */
void TimeUpdate(Matrix *P, Matrix *Phi, Matrix *Qdt=0);

#endif //C___TIMEUPDATE_H
