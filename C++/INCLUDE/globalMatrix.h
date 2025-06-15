#ifndef GLOBALMATRIX_H
#define GLOBALMATRIX_H

#include "matrix.h"
#include <Mjday.h>
#include <SAT_Const.h>
#include <string.h>

/**
 * @brief Matriz de coeficientes Cnm del campo gravitacional armónico.
 *
 * Dimensiones: 181 x 181
 */
extern Matrix Cnm;

/**
 * @brief Matriz de coeficientes Snm del campo gravitacional armónico.
 *
 * Dimensiones: 181 x 181
 */
extern Matrix Snm;

/**
 * @brief Matriz de observaciones del satélite GEOS3.
 *
 * Columnas:
 * - [1] Tiempo (MJD)
 * - [2] Azimut [rad]
 * - [3] Elevación [rad]
 * - [4] Distancia [m]
 *
 * Dimensiones: 46 x 4
 */
extern Matrix obs;

extern Matrix eopdata;


/**
 * @brief Carga el modelo gravitacional GGM03S desde un archivo de datos.
 *
 * Lee el archivo `../DATA/GGM03S.txt` y carga los coeficientes armónicos
 * del campo gravitacional en las matrices globales `Cnm` y `Snm`.
 *
 * @note Se espera que el archivo contenga 181 filas y columnas según el modelo.
 * @throws Sale del programa si el archivo no puede abrirse.
 */
void GGM03S();

/**
 * @brief Carga las observaciones del satélite GEOS3.
 *
 * Lee el archivo `../DATA/GEOS3.txt` y almacena los datos de tiempo,
 * azimut, elevación y distancia en la matriz global `obs`.
 *
 * @throws Sale del programa si el archivo no puede abrirse.
 */
void GEOS3();

void EOPDATA();

#endif // GLOBALMATRIX_H
