
#include "../INCLUDE/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

Matrix::Matrix(int fil, int col): fil(fil),col(col) { //fil(fil) da valores a la funcion miembro fil, equivalente a decir this->fil=fil
    initMatrix();
}
Matrix::Matrix(int fil, int col, const double *v, int n): fil(fil),col(col) {
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
    if (this == &matrix2) return *this;

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

Matrix Matrix::operator*(double hi) const {
    Matrix result(fil, col);
    for (int i = 0; i < fil; i++) {
        for (int j = 0; j < col; j++) {
            result.matrix[i][j] = matrix[i][j] * hi;
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
Matrix Matrix::transpose() const{
    Matrix resultado(this->col, this->fil);
    for (int i=1; i<=this->fil;i++) {
        for (int j=1; j<=this->col;j++) {
            resultado(j,i)=(*this)(i,j);
        }
    }
    return resultado;
}

Matrix Matrix::identity() {
    Matrix resultado(fil, col);
    for (int i=1; i<=fil; i++) {
        for (int j=1; j<=col; j++) {
            if (i==j) {
                resultado(i,j) = 1;
            }
            else {
                resultado(i,j) = 0;
            }
        }
    }
    return resultado;
}

Matrix Matrix::inverse() const {
    int n = fil;
    Matrix A(fil, col);
    // Copiar los valores de la matriz original a A usando un bucle
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            A(i, j) = (*this)(i, j);
        }
    }
    Matrix I(n, n);
    I = I.identity();
    for (int i = 1; i <= n; i++) {
        double maxVal = fabs(A(i, i));
        int pivoteFila = i;
        for (int k = i + 1; k <= n; k++) {
            if (fabs(A(k, i)) > maxVal) {
                maxVal = fabs(A(k, i));
                pivoteFila = k;
            }
        }
        if (pivoteFila != i) {
            for (int j = 1; j <= n; j++) {
                swap(A(i, j), A(pivoteFila, j));
                swap(I(i, j), I(pivoteFila, j));
            }
        }
        double diagonal = A(i, i);
        for (int j = 1; j <= n; j++) {
            A(i, j) = A(i, j) / diagonal;
            I(i, j) = I(i, j) / diagonal;
        }
        for (int k = 1; k <= n; ++k) {
            if (k == i) continue;
            double factor = A(k, i);
            for (int j = 1; j <= n; ++j) {
                A(k, j) -= factor * A(i, j);
                I(k, j) -= factor * I(i, j);
            }
        }
    }
    return I;
}
int Matrix::getFil() const {
    return fil;
}
int Matrix::getCol() const {
    return col;
}



Matrix Matrix::LoadFromFile(string filename){
    string path = "../DATA/" + filename;
    ifstream file(path);
    if(!file.is_open()){
        throw runtime_error("No existe el archivo");
    }
    string s;
    int rows = 0;
    int cols = 0;

    while (std::getline(file, s)) {
        std::istringstream stream(s);
        double value;
        int currentCols = 0;

        // Contamos el número de columnas en la fila
        while (stream >> value) {
            currentCols++;
        }

        if (cols == 0) {
            cols = currentCols;  // Asumimos que todas las filas tienen el mismo número de columnas
        } else if (currentCols != cols) {
            throw std::runtime_error("Las filas tienen diferentes números de columnas.");
        }

        rows++;
    }

    file.clear();
    file.seekg(0, std::ios::beg);

    Matrix mat(rows, cols);

    int i = 1;
    while (std::getline(file, s)) {
        std::istringstream stream(s);
        double value;
        int j = 1;

        while (stream >> value) {
            mat(i, j) = value;
            j++;
        }
        i++;
    }

    file.close();
    return mat;
}

double Matrix::norm() const {
    double suma = 0;
    for (int i = 1; i <= fil; i++) {
        for (int j = 1; j <= fil; j++) {
            suma += (*this)(i, j) * (*this)(i, j);
        }
    }
    return sqrt(suma);
}
double Matrix::normVec() const {
    double suma = 0;
    for (int i = 1; i <= fil; i++) {
        for (int j = 1; j <= col; j++) {
            suma += (*this)(i, j) * (*this)(i, j);
        }
    }
    return sqrt(suma);
}






