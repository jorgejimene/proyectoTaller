//
// Created by jojimene on 23/04/2025.
//
/**
 * @file  SAT_Const.h
 * @brief Constantes físicas y astronómicas usadas en el proyecto
 **********************/

#ifndef C___SAT_CONST_H
#define C___SAT_CONST_H
#include <cmath>
/**
 * @struct Constants
 * @brief Contiene constantes matemáticas, físicas y astronómicas para cálculos satelitales
 */
struct Constants {
    const double pi2;          ///< 2 * pi
    const double Rad;          ///< Factor de conversión grados a radianes
    const double Deg;          ///< Factor de conversión radianes a grados
    const double Arcs;         ///< Factor de conversión de segundos de arco a radianes

    const double MJD_J2000;    ///< Fecha Juliana Modificada del epoch J2000
    const double T_B1950;      ///< Fecha de referencia B1950
    const double c_light;      ///< Velocidad de la luz (m/s)
    const double AU;           ///< Unidad Astronómica (m)

    const double R_Earth;      ///< Radio ecuatorial de la Tierra (m)
    const double f_Earth;      ///< Aplanamiento de la Tierra
    const double R_Sun;        ///< Radio del Sol (m)
    const double R_Moon;       ///< Radio de la Luna (m)

    const double omega_Earth;  ///< Velocidad angular de rotación de la Tierra (rad/s)

    const double GM_Earth;     ///< Producto gravitacional estándar de la Tierra (m³/s²)
    const double GM_Sun;       ///< Producto gravitacional estándar del Sol (m³/s²)
    const double GM_Moon;      ///< Producto gravitacional estándar de la Luna (m³/s²)
    const double GM_Mercury;   ///< Producto gravitacional estándar de Mercurio (m³/s²)
    const double GM_Venus;     ///< Producto gravitacional estándar de Venus (m³/s²)
    const double GM_Mars;      ///< Producto gravitacional estándar de Marte (m³/s²)
    const double GM_Jupiter;   ///< Producto gravitacional estándar de Júpiter (m³/s²)
    const double GM_Saturn;    ///< Producto gravitacional estándar de Saturno (m³/s²)
    const double GM_Uranus;    ///< Producto gravitacional estándar de Urano (m³/s²)
    const double GM_Neptune;   ///< Producto gravitacional estándar de Neptuno (m³/s²)
    const double GM_Pluto;     ///< Producto gravitacional estándar de Plutón (m³/s²)

    const double P_Sol;        ///< Presión de radiación solar a 1 AU (N/m²)
};

// Declaración de una instancia global de la estructura Constants
extern const Constants consts;

#endif //C___SAT_CONST_H
