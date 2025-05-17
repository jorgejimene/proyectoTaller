//
// Created by jojimene on 07/05/2025.
//

#ifndef C___JPL_EPH_DE430_H
#define C___JPL_EPH_DE430_H
void JPL_Eph_DE430(double Mjd_TDB,double r_Mercury[3],double r_Venus[3],double r_Earth[3],double r_Mars[3],double r_Jupiter[3],double r_Saturn[3],double r_Uranus[3], double r_Neptune[3],double r_Pluto[3],double r_Moon[3],double r_Sun[3]);
void calcular_posicion(const double* PCtemp, int base_idx, int grado, int j, double Mjd_TDB, double Mjd0, double r_out[3]);
#endif //C___JPL_EPH_DE430_H
