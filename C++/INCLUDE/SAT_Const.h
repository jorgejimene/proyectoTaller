//
// Created by jojimene on 23/04/2025.
//

#ifndef C___SAT_CONST_H
#define C___SAT_CONST_H

struct Constants {
    const double pi2;
    const double Rad;
    const double Deg;
    const double Arcs;

    const double MJD_J2000;
    const double T_B1950;
    const double c_light;
    const double AU;

    const double R_Earth;
    const double f_Earth;
    const double R_Sun;
    const double R_Moon;

    const double omega_Earth;

    const double GM_Earth;
    const double GM_Sun;
    const double GM_Moon;
    const double GM_Mercury;
    const double GM_Venus;
    const double GM_Mars;
    const double GM_Jupiter;
    const double GM_Saturn;
    const double GM_Uranus;
    const double GM_Neptune;
    const double GM_Pluto;

    const double P_Sol;
};

// Declaración de una instancia global de la estructura Constants
extern const Constants consts;

#endif //C___SAT_CONST_H
