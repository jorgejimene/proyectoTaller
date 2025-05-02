#include "../INCLUDE/MeasUpdate.h"
#include "../INCLUDE/Matrix.h"

#include "../INCLUDE/MeasUpdate.h"
#include "../INCLUDE/Matrix.h"

void MeasUpdate(Matrix &x, Matrix z, const Matrix &g, const Matrix &s,
                Matrix G, Matrix &P, int n, Matrix &K) {
    int m = z.getFil();

    // 1. Create inverse weight matrix
    Matrix Inv_W(m,m);
    for (int i=1; i<m; i++) {
        Inv_W(i,i) = s(i,0)*s(i,0);  // Note: changed to 0-based indexing
    }

    // 2. Calculate Kalman gain
    Matrix GTras = G.transpose();
    Matrix temp = G*P*GTras;
    Matrix inversa = (Inv_W + temp).inverse();
    K = P*GTras*inversa;

    // 3. State update
    x = x + K*(z - g);

    // 4. Covariance update
    Matrix I = Matrix(n,n).identity();
    P = (I - K*G)*P;
}