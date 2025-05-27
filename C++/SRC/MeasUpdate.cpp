#include "../INCLUDE/MeasUpdate.h"

void MeasUpdate(Matrix &x, Matrix &z, Matrix &g, Matrix &s,
                Matrix &G, Matrix &P, int n, Matrix &K) {
    int m = z.getFil();


    Matrix Inv_W = Matrix(m,m);
    for (int i = 1; i <= m; i++) {
        Inv_W(i, i) = s(i, 1) * s(i, 1);
    }

    Matrix GTras = G.transpose();
    Matrix temp = G * P * GTras;
    Matrix suma = Inv_W + temp;
    Matrix inversa = suma.inverse();
    K = P * GTras * inversa;
    Matrix x_new = x + K*(z-g);
    x = x_new;

    Matrix I(n, n);
    I = I.identity();
    P = (I - K * G) * P;
}
