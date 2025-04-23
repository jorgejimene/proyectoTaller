
#include "../INCLUDE/Matrix.h"
#include <iostream>
#include <iomanip>

using namespace std;

Matrix::Matrix(int fil, int col): fil(fil),col(col) { //fil(fil) da valores a la funcion miembro fil, equivalente a decir this->fil=fil
    initMatrix();
}
Matrix::Matrix(int fil, int col, double *v, int n): fil(fil),col(col) {
    initMatrix();

    int k = 0;

    for (int i=0; i<fil; i++){
        for(int j=0;j<col;j++){
            if(k<n){
                matrix[i][j]=v[k++];
            }else{
                matrix[i][j]=0;
            }
        }
    }
}


Matrix::Matrix(const Matrix &m) {
    *this = m;
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
Matrix::~Matrix() {
    for (int i=0; i<fil; i++){
        delete[] matrix[i];
    }
    delete[] matrix;
}
Matrix &Matrix::operator=(const Matrix &matrix2) {
    for(int i=0; i< fil; i++){
        for(int j=0; j<col;j++){
            this->matrix[i][j] = matrix2.matrix[i][j];
        }
    }
    return *this;
}
Matrix Matrix::operator+(const Matrix &matrix2) {
    Matrix result(fil, col);

    for(int i=0; i < fil; i++){
        for(int j=0;j<col;j++){
            result.matrix[i][j]=matrix[i][j] + matrix2.matrix[i][j];
        }
    }
    return result;
}
Matrix Matrix::operator-(const Matrix &matrix2) {
    Matrix result(fil, col);

    for(int i=0; i < fil; i++){
        for(int j=0;j<col;j++){
            result.matrix[i][j]=matrix[i][j] - matrix2.matrix[i][j];
        }
    }
    return result;
}
Matrix Matrix::operator*(const Matrix &matrix2) {
    Matrix result(fil, col);

    for(int i=0; i < this->fil; i++){
        for(int j=0; j<matrix2.col; j++){
            result.matrix[i][j]=0;
            for(int k=0; k<this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
    return result;
}
double &Matrix::operator()(const int i, const int j) const {
    return matrix[i-1][j-1];
}
void Matrix::print() {
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            cout << fixed << setprecision(14) << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
void Matrix::printMatrixValues(float** arr, int n, int m){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout<<arr[i][j]<<"\t";
        }
        cout<<endl;
    }
    return;
}
void Matrix::printInverseMatrix(float** arr, int n, int m){
    for (int i = 0; i < n; i++) {
        for (int j = n; j < m; j++) {
            printf("%.3f\t", arr[i][j]);
        }
        cout<<endl;
    }
    return;
}
void Matrix::findInvMatGaussJordan(float** mat, int order){
    float temp;
    printf("The inverse of matrix : A = \n");
    printMatrixValues(mat, order, order);
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < 2 * order; j++) {
            if (j == (i + order))
                mat[i][j] = 1;
        }
    }
    for (int i = order - 1; i > 0; i--) {
        if (mat[i - 1][0] < mat[i][0]) {
            float* temp = mat[i];
            mat[i] = mat[i - 1];
            mat[i - 1] = temp;
        }
    }
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            if (j != i) {
                temp = mat[j][i] / mat[i][i];
                for (int k = 0; k < 2 * order; k++) {
                    mat[j][k] -= mat[i][k] * temp;
                }
            }
        }
    }
    for (int i = 0; i < order; i++) {
        temp = mat[i][i];
        for (int j = 0; j < 2 * order; j++) {
            mat[i][j] = mat[i][j] / temp;
        }
    }
    cout<<"A' =\n";
    printInverseMatrix(mat, order, 2 * order);
    return;
}
