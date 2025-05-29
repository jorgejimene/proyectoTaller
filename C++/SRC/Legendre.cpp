#include <cmath>
#include "../INCLUDE/Matrix.h"

using namespace std;

void Legendre(int n, int m, double fi, Matrix &pnm, Matrix &dpnm) {

    pnm(1, 1) = 1.0;
    dpnm(1, 1) = 0.0;

    if (n >= 1 && m >= 1) {
        pnm(2, 2) = sqrt(3.0) * cos(fi);
        dpnm(2, 2) = -sqrt(3.0) * sin(fi);
    }

    // Diagonal coefficients
    for (int i = 2; i <= n; ++i) {
        pnm(i + 1, i + 1) = sqrt((2.0 * i + 1.0) / (2.0 * i)) * cos(fi) * pnm(i, i);
    }

    for (int i = 2; i <= n; ++i) {
        dpnm(i + 1, i + 1) = sqrt((2.0 * i + 1.0) / (2.0 * i)) *
                            (cos(fi) * dpnm(i, i) - sin(fi) * pnm(i, i));
    }

    // Horizontal first step coefficients
    for (int i = 1; i <= n; ++i) {
        pnm(i + 1, i) = sqrt(2.0 * i + 1.0) * sin(fi) * pnm(i, i);
    }

    for (int i = 1; i <= n; ++i) {
        dpnm(i + 1, i) = sqrt(2.0 * i + 1.0) *
                         (cos(fi) * pnm(i, i) + sin(fi) * dpnm(i, i));
    }

    // Horizontal second step coefficients
    int j = 0;
    int k = 2;
    double a,b,c,d;
    while (true) {
        for (int i = k; i <= n; ++i) {
            a = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j)));
            b = sqrt(2.0 * i - 1.0) * sin(fi) * pnm(i, j + 1);
            c = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * pnm(i - 1, j + 1);
            pnm(i + 1, j + 1) = a * (b - c);
        }
        ++j;
        ++k;
        if (j > m) break;
    }

    j = 0;
    k = 2;

    while (true) {
        for (int i = k; i <= n; ++i) {
            a = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j)));
            b = sqrt(2.0 * i - 1.0) * sin(fi) * dpnm(i, j + 1);
            c = sqrt(2.0 * i - 1.0) * cos(fi) * pnm(i, j + 1);
            d = sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * dpnm(i - 1, j + 1);
            dpnm(i + 1, j + 1) = a * (b + c - d);
        }
        ++j;
        ++k;
        if (j > m) break;
    }
}
