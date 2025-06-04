//
// Created by jojimene on 03/04/2025.
//
/**
 * @file  Mjday.h
 * @brief Función Mjday usada en el proyecto, transforma fechas a juliano
 **********************/

#ifndef C___MJDAY_H
#define C___MJDAY_H
#include <cmath>
/**
 * @details Convierte una fecha y hora en formato calendario (año, mes, día, hora, minuto, segundo)
 * a su equivalente en MJD
 * @param yr año
 * @param mon mes
 * @param day dia
 * @param hr hora
 * @param min minuto
 * @param sec segundos
 * @return Valor de la fecha juliana
 */

double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0.0);

#endif //C___MJDAY_H
