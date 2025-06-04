/**
* @file Cheb3D.h
 * @brief Aproximación de vectores tridimensionales usando polinomios de Chebyshev.
 */

#ifndef CHEB3D_H
#define CHEB3D_H

/**
 * @brief Calcula la aproximación de Chebyshev de un vector 3D en un instante dado.
 *
 * @param t        Tiempo de evaluación.
 * @param N        Número de coeficientes del polinomio.
 * @param Ta       Inicio del intervalo de tiempo.
 * @param Tb       Fin del intervalo de tiempo.
 * @param Cx       Coeficientes de Chebyshev para la coordenada X.
 * @param Cy       Coeficientes de Chebyshev para la coordenada Y.
 * @param Cz       Coeficientes de Chebyshev para la coordenada Z.
 * @param ChebApp  Vector de salida con la aproximación (dimensión 3).
 *
 * @throws std::out_of_range si el tiempo t está fuera del intervalo [Ta, Tb].
 */
void Cheb3D(double t, int N, double Ta, double Tb,
            const double Cx[], const double Cy[], const double Cz[],
            double ChebApp[3]);

#endif // CHEB3D_H
