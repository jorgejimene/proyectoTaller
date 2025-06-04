#ifndef GAST_H
#define GAST_H
#include "gmst.h"
#include "EqnEquinox.h"
#include "cmath"
/**
 * @brief Calcula el Tiempo Sidéreo Aparente de Greenwich (GAST).
 *
 *
 * @param Mjd_UT1 Fecha juliana modificada (UT1).
 * @return Tiempo sidéreo aparente de Greenwich en radianes.
 */
double gast(double Mjd_UT1);

#endif // GAST_H
