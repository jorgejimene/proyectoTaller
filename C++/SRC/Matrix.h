//
// Created by jojimene on 20/03/2025.
//

#ifndef C___MATRIX_H
#define C___MATRIX_H


class Matrix {
public:
    Matrix(int fil, int col); //reserva memoria
    Matrix(int fil, int col, double v[], int n); //reserva memoria y permite introducir datos almacenados en un vector, siendo n la dimensión del vector
    Matrix(const Matrix& m);
    ~Matrix(); //destructor

    //sobrecarga de los operadores
    //el & en este caso significa "referencia", mecanismo de c++ que permite que a una funcion le podamos asignar algo
    Matrix& operator=(const Matrix& matrix2);
    Matrix operator+(const Matrix& matrix2);
    Matrix operator-(const Matrix& matrix2);
    Matrix operator*(const Matrix& matrix2);
    //para acceder a los valores de la matriz, por ejemplo a(i,j)
    double& operator()(const int i, const int j)
        const;

    void print();

private:
    void initMatrix();

private:
    int fil;
    int col;
    double **matrix;
};


#endif //C___MATRIX_H
