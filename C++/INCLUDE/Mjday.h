//
// Created by jojimene on 03/04/2025.
//
/**
 * @file  Mjday.h
 * @author Jorge
 * @brief Clase Mjday usada en el proyecto, transforma fechas a juliano
 * @date 9-4-2025
 **********************/

#ifndef C___MJDAY_H
#define C___MJDAY_H
/**
 * @details Conversor de fechas a sus valores julianos
 * @param yr año
 * @param mon mes
 * @param day dia
 * @param hr hora
 * @param min minuto
 * @param sec segundos
 * @return Valor de la fecha en juliano
 */

double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0.0);

#endif //C___MJDAY_H
