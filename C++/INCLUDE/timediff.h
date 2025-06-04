//
// Created by jojimene on 07/05/2025.
//
/**
 * @file  timediff.h
 * @brief Función timediff para calcular diferencias de tiempo entre UT1, UTC, TAI, GPS y TT
 **********************/

#ifndef C___TIMEDIFF_H
#define C___TIMEDIFF_H

/**
 * @brief Calcula diferencias de tiempo entre distintas escalas temporales
 *
 * @param UT1_UTC Diferencia UT1-UTC [s]
 * @param TAI_UTC Diferencia TAI-UTC [s]
 * @param UT1_TAI Referencia para devolver UT1-TAI [s]
 * @param UTC_GPS Referencia para devolver UTC-GPS [s]
 * @param UT1_GPS Referencia para devolver UT1-GPS [s]
 * @param TT_UTC Referencia para devolver TT-UTC [s]
 * @param GPS_UTC Referencia para devolver GPS-UTC [s]
 */
void timediff(double UT1_UTC, double TAI_UTC, double &UT1_TAI, double &UTC_GPS, double &UT1_GPS, double &TT_UTC, double &GPS_UTC);

#endif //C___TIMEDIFF_H
