#include "../INCLUDE/MeasUpdate.h"
#include "../INCLUDE/Matrix.h"

void MeasUpdate(Matrix &x, Matrix z, Matrix g,Matrix s, Matrix G, Matrix &P, int n, Matrix &K) {
    int m = z.getFil();

    Matrix Inv_W(m,m);
    Matrix A(n,n);

    for (int i=1; i<=m;i++) {
        Inv_W(i,i) = s(i,1)*s(i,1);    // Inverse weight (measurement covariance)
    }

    // Kalman gain
    Matrix GTras = G.transpose();
    GTras.print();
    Matrix inversa= (Inv_W+G*P*GTras).inverse();

    K = P*GTras*inversa;

    // State update
    x = x + K*(z-g);

    // Covariance update
    A = A.identity();
    P = (A-K*G)*P;
}
