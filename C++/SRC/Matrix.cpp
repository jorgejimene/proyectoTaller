//
// Created by jojimene on 20/03/2025.
//

#include "../INCLUDE/Matrix.h"
#include <iostream>

Matrix::Matrix(int fil, int col): fil(fil),col(col) {
    initMatrix();
}

Matrix::Matrix(int fil, int col, double *v, int n): fil(fil),col(col) {
    initMatrix();
}



Matrix::Matrix(const Matrix &m) {

}

void Matrix::initMatrix() {
    //reserva de memoria de las columnas
    matrix = new double*[fil];
    for(int i=0;i<fil;i++){
        matrix[i]=new double[col];
    }
    //añadimos valores a las filas
    for (int i =0;i<fil;i++){
        for(int j=0;j<col;j++){
            matrix[i][j]=0.0;
        }
    }
}