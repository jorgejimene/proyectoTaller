//
// Created by jojimene on 20/03/2025.
//

/**
 * @file  Matrix.h
 * @author Jorge
 * @brief Clase matriz usada en el proyecto
 * @date 26-03-2025
 **********************/

#ifndef C___MATRIX_H
#define C___MATRIX_H


class Matrix {
public:
    int getFil() const;
    int getCol() const;

    /**
     * @details Constructor básico de la clase Matriz, crea una matriz rellenada con ceros
     * @param fil Número de filas de la matriz
     * @param col Número de columnas de la matriz
     * @return Objeto matriz con fil filas y col columnas
     */

    Matrix(int fil, int col); //reserva memoria
    /**
     * @details Constructor de la clase matriz que construye una matriz con los valores de v
     * fil numero de filas y col numero de columnas
     * @param fil Numero de filas de la matriz
     * @param col Numero de columnas de la matriz
     * @param v Vector con los valores a introducir en la matriz
     * @param n Tamaño del vector v
     * @return Objeto matriz con fil filas, col columnas y los valores de v
     */
    Matrix(int fil, int col, double v[], int n); //reserva memoria y permite introducir datos almacenados en un vector, siendo n la dimensión del vector
    /**
     * @details Función que hace que *this apunte a la matriz m
     * @param m Matriz a la que queremos apuntar
     */
    Matrix(const Matrix& m);
    /**
     * @brief Destructor de la clase matriz
     */
    ~Matrix(); //destructor

    //sobrecarga de los operadores
    //el & en este caso significa "referencia", mecanismo de c++ que permite que a una funcion le podamos asignar algo
    /**
     * @details Sobrecarga del operador =
     * @param matrix2  Matriz a la que queremos igualar
     * @return Objeto matriz igualada a los valores de matrix2
     */
    Matrix& operator=(const Matrix& matrix2);
    /**
     * @details Sobrecarga del operador +
     * @param matrix2 Matriz con la que queremos sumar
     * @return Matriz resultante de la suma
     */
    Matrix operator+(const Matrix& matrix2);
    /**
     * @details Sobrecarga del operador -
     * @param matrix2 Matriz con la que queremos restar
     * @return Matriz resultante de la resta
     */
    Matrix operator-(const Matrix& matrix2);
    /**
     * @details Sobrecarga del operador *
     * @param matrix2 Matriz con la que queremos multiplicar
     * @return Matriz resultante del producto
     */
    Matrix operator*(const Matrix& matrix2);
    /**
     * @details Sobrecarga del operador ()
     * @param i Valor de la fila donde queremos buscar
     * @param j Valor de la columna donde queremos buscar
     * @return Devuelve el valor en la posición (i,j) de la matiz
     */
    //para acceder a los valores de la matriz, por ejemplo a(i,j)
    double& operator()(const int i, const int j) const; //este const se ha puesto para "proteger" y que i y j no se cambien
    /**
     * @details Muestra por pantalla la matriz
     */
    void print();
    Matrix transpose() const;
    Matrix identity();
    Matrix inverse() const;


private:
    int fil;
    int col;
    /**
     * @details Función de inicialización de la matriz
     */
    void initMatrix(); //funcion miembro que solo voy a poder usar yo como programador, ahí se reservará la memoria

    double **matrix;
    void findInvMatGaussJordan(float **mat, int order);


    void printMatrixValues(float **arr, int n, int m);
    void printInverseMatrix(float ** arr, int n, int m);

};


#endif //C___MATRIX_H
