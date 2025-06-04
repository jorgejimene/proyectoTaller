
#include "../INCLUDE/SAT_Const.h"

const Constants consts = {
        // Mathematical constants
        2 * M_PI,              // pi2
        M_PI / 180,            // Rad
        180 / M_PI,            // Deg
        3600 * 180 / M_PI,     // Arcs

        // General
        51544.5,               // MJD_J2000
        -0.500002108,          // T_B1950
        299792458.0,           // c_light
        149597870700.0,        // AU

        // Physical parameters of the Earth, Sun, and Moon
        6378.1363e3,           // R_Earth
        1 / 298.257223563,     // f_Earth
        696000e3,              // R_Sun
        1738e3,                // R_Moon

        // Earth rotation (omega_Earth)
        15.04106717866910 / 3600 * (M_PI / 180), // omega_Earth

        // Gravitational coefficients
        398600.435436e9,       // GM_Earth
        132712440041.939400e9, // GM_Sun
        398600.435436e9 / 81.30056907419062, // GM_Moon
        22031.780000e9,        // GM_Mercury
        324858.592000e9,       // GM_Venus
        42828.375214e9,        // GM_Mars
        126712764.800000e9,    // GM_Jupiter
        37940585.200000e9,     // GM_Saturn
        5794548.600000e9,      // GM_Uranus
        6836527.100580e9,      // GM_Neptune
        977.0000000000009e9,   // GM_Pluto

        // Solar radiation pressure at 1 AU
        1367 / 299792458.0     // P_Sol
};