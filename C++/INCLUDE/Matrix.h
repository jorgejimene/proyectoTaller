//
// Created by jojimene on 20/03/2025.
//

/**
 * @file  Matrix.h
 * @brief Clase matriz usada en el proyecto
 **********************/

#ifndef C___MATRIX_H
#define C___MATRIX_H

#include <string>

class Matrix {
public:
    /**
     * @brief Obtiene el número de filas de la matriz.
     * @return Número de filas.
     */
    int getFil() const;

    /**
     * @brief Obtiene el número de columnas de la matriz.
     * @return Número de columnas.
     */
    int getCol() const;

    static Matrix & transposeDif(Matrix &m);


    /**
     * @brief Constructor por defecto. Inicializa una matriz de ceros.
     * @param fil Número de filas.
     * @param col Número de columnas.
     */
    Matrix(int fil, int col);

    /**
     * @brief Constructor con datos.
     * @param fil Número de filas.
     * @param col Número de columnas.
     * @param v Vector de datos.
     * @param n Tamaño del vector.
     */
    Matrix(int fil, int col, double v[], int n);

    /**
     * @brief Constructor de copia.
     * @param m Matriz a copiar.
     */
    Matrix(const Matrix& m);

    /**
     * @brief Destructor.
     */
    ~Matrix();

    /**
     * @brief Sobrecarga del operador de asignación.
     * @param matrix2 Matriz a asignar.
     * @return Referencia a esta matriz.
     */
    Matrix& operator=(const Matrix& matrix2);

    /**
     * @brief Sobrecarga del operador suma.
     * @param matrix2 Matriz a sumar.
     * @return Matriz resultante.
     */
    Matrix operator+(const Matrix& matrix2);

    /**
     * @brief Sobrecarga del operador resta.
     * @param matrix2 Matriz a restar.
     * @return Matriz resultante.
     */
    Matrix operator-(const Matrix& matrix2);

    /**
     * @brief Sobrecarga del operador multiplicación.
     * @param matrix2 Matriz a multiplicar.
     * @return Matriz resultante.
     */
    Matrix operator*(const Matrix& matrix2);
    /**
     * @brief Sobrecarga del operador multiplicación con double.
     * @param hi Double a mutliplicar.
     * @return Matriz resultante.
     */
    Matrix operator*(double hi) const;

    /**
     * @brief Accede al valor de una posición de la matriz.
     * @param i Fila.
     * @param j Columna.
     * @return Referencia al valor (i,j).
     */
    double& operator()(const int i, const int j) const;

    /**
     * @brief Imprime la matriz por pantalla.
     */
    void print();

    /**
     * @brief Devuelve la matriz transpuesta.
     * @return Matriz transpuesta.
     */
    Matrix transpose() const;

    /**
     * @brief Devuelve una matriz identidad del mismo tamaño.
     * @return Matriz identidad.
     */
    Matrix identity();

    /**
     * @brief Devuelve la matriz inversa.
     * @return Matriz inversa.
     */
    Matrix inverse() const;

    /**
     * @brief Carga una matriz desde archivo.
     * @param filename Nombre del archivo.
     * @return Matriz cargada.
     */
    static Matrix LoadFromFile(std::string filename);

    /**
     * @brief Calcula la norma (matricial) de la matriz.
     * @return Norma.
     */
    double norm() const;

    /**
     * @brief Calcula la norma del vector (si aplica).
     * @return Norma del vector.
     */
    double normVec() const;

private:
    int fil;  ///< Número de filas.
    int col;  ///< Número de columnas.

    /**
     * @brief Inicializa la matriz con memoria dinámica.
     */
    void initMatrix();

    double **matrix;  ///< Puntero doble para la matriz.

    /**
     * @brief Imprime los valores de una matriz tipo float.
     * @param arr Puntero a la matriz.
     * @param n Filas.
     * @param m Columnas.
     */
    void printMatrixValues(float **arr, int n, int m);

    /**
     * @brief Imprime la matriz inversa desde una matriz aumentada.
     * @param arr Matriz aumentada.
     * @param n Número de filas.
     * @param m Número de columnas.
     */
    void printInverseMatrix(float ** arr, int n, int m);


};

#endif //C___MATRIX_H
