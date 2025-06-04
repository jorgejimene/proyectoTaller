//
// Created by jojimene on 07/05/2025.
//
/**
 * @file  Position.h
 * @brief Convierte coordenadas geodésicas (longitud, latitud, altura) a vector de posición en metros
 **********************/

#ifndef C___POSITION_H
#define C___POSITION_H
#include <cmath>
#include "SAT_Const.h"
/**
 * @details Calcula el vector de posición (r) en coordenadas cartesianas a partir
 *          de coordenadas geodésicas: longitud (rad), latitud (rad) y altura (m).
 * @param lon Longitud en radianes
 * @param lat Latitud en radianes
 * @param h Altura sobre el elipsoide en metros
 * @param r Puntero a arreglo de 3 elementos para almacenar el vector posición [m]
 */

void Position(double lon, double lat, double h, double* r);

#endif //C___POSITION_H
