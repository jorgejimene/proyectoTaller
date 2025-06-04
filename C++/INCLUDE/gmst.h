#ifndef GMST_H
#define GMST_H
#include <math.h>
#include "Frac.h"
/**
 * @brief Calcula el Tiempo Sideral Medio de Greenwich (GMST).
 *
 *
 * @param Mjd_UT1 Fecha Juliana Modificada (UT1)
 * @return Tiempo Sideral Medio de Greenwich en radianes [0, 2π]
 */
double gmst(double Mjd_UT1);

#endif // GMST_H
