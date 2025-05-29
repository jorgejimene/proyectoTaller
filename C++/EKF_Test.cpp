/*
 * Archivo de tests
 *
 */

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "Accel.h"
#include "AccelHarmonic.h"
#include "AuxParam.h"
#include "Cheb3D.h"
#include "globalMatrix.h"
#include "G_AccelHarmonic.h"
#include "VarEqn.h"
#include "./INCLUDE/Matrix.h"
#include "./INCLUDE/Mjday.h"
#include "./INCLUDE/R_x.h"
#include "./INCLUDE/R_y.h"
#include "./INCLUDE/R_z.h"
#include "./INCLUDE/Legendre.h"
#include "./INCLUDE/SAT_Const.h"
#include "./INCLUDE/NutAngles.h"
#include "./INCLUDE/Unit.h"
#include "./INCLUDE/AzElPa.h"
#include "./INCLUDE/Frac.h"
#include "./INCLUDE/AccelPointMass.h"
#include "./INCLUDE/IERS.h"
#include "./INCLUDE/MeanObliquity.h"
#include "./INCLUDE/TimeUpdate.h"
#include "./INCLUDE/Mjday_TBD.h"
#include "./INCLUDE/MeasUpdate.h"
#include "./INCLUDE/Sign_.h"
#include "./INCLUDE/NutMatrix.h"
#include "./INCLUDE/PoleMatrix.h"
#include "./INCLUDE/PrecMatrix.h"
#include "./INCLUDE/timediff.h"
#include "./INCLUDE/EccAnom.h"
#include "./INCLUDE/elements.h"
#include "./INCLUDE/EqnEquinox.h"
#include "./INCLUDE/gmst.h"
#include "./INCLUDE/gast.h"
#include "./INCLUDE/GHAMatrix.h"
#include "./INCLUDE/JPL_Eph_DE430.h"
#include "./INCLUDE/LTC.h"


/*
#include "./SRC/Matrix.cpp"
#include "./SRC/Legendre.cpp"
#include "./SRC/SAT_Const.cpp"
#include "./SRC/NutAngles.cpp"
#include "./SRC/Unit.cpp"
#include "./SRC/AzElPa.cpp"
#include "./SRC/Frac.cpp"
#include "./SRC/Mjday.cpp"
#include "./SRC/R_x.cpp"
#include "./SRC/R_y.cpp"
#include "./SRC/R_z.cpp"
#include "./SRC/AccelPointMass.cpp"
#include "./SRC/IERS.cpp"
 #include "./SRC/MeanObliquity.cpp"
*/

int tests_run = 0;

#define TOL_ 10e-14
#define TOL2_ 10e-10 //para numeros mas pequeños
#define TOL3_ 10e-9
#define TOL4_ 10e-3

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;

int testConstructorVacio(){

    Matrix m1(2,2);
    _assert(m1(1,1)==0 && m1(1,2)==0 && m1(2,1)==0 && m1(2,2)==0);
    return 0;
}
int testConstructorVector(){
    double v[] = {1.0, 2.0, 3.0, 4.0};
    Matrix m1(2,2,v,4);
    _assert(m1(1,1)==1 && m1(1,2)==2 && m1(2,1)==3 && m1(2,2)==4);
    return 0;
}
int testSuma(){
    double v[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {5.0, 6.0, 7.0, 8.0};

    double resultado[] = {6,8,10,12};

    Matrix m1 (2,2,v,4);
    Matrix m2(2,2,v2,4);
    Matrix m3 = m1 + m2;

    Matrix m4 (2,2,resultado,4);
    _assert(m3(1,1)==m4(1,1) && m3(1,2)==m4(1,2)
            && m3(2,1)==m4(2,1) && m3(2,2)==m4(2,2));
    return 0;
}
int testResta(){
    double v[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {5.0, 6.0, 7.0, 8.0};

    double resultado[] = {4,4,4,4};

    Matrix m1 (2,2,v,4);
    Matrix m2(2,2,v2,4);
    Matrix m3 = m2 - m1;

    Matrix m4 (2,2,resultado,4);
    _assert(m3(1,1)==m4(1,1) && m3(1,2)==m4(1,2)
            && m3(2,1)==m4(2,1) && m3(2,2)==m4(2,2));
    return 0;
}
int testProducto(){
    double v[] = {1.0, 2.0, 3.0, 4.0};
    double v2[] = {5.0, 6.0, 7.0, 8.0};

    double resultado[] = {19,22,43,50};

    Matrix m1 (2,2,v,4);
    Matrix m2(2,2,v2,4);
    Matrix m3 = m1 * m2;

    Matrix m4 (2,2,resultado,4);
    _assert(m3(1,1)==m4(1,1) && m3(1,2)==m4(1,2)
            && m3(2,1)==m4(2,1) && m3(2,2)==m4(2,2));
    return 0;
}
int testTraspuesta() {
    double valores[] = {
        1, 2, 3,
        4, 5, 6,
    };
    Matrix m(2,3,valores,6);
    Matrix traspM = m.transpose();
    traspM.print();

    _assert(traspM.getCol() == 2);
    _assert(traspM.getFil() == 3);
    _assert(traspM(1, 1) == 1);
    _assert(traspM(2, 1) == 2);
    _assert(traspM(3, 1) == 3);
    _assert(traspM(1, 2) == 4);
    _assert(traspM(2, 2) == 5);
    _assert(traspM(3, 2) == 6);

    return 0;
}
int Mjday_01(){

    _assert(fabs(Mjday(2025, 4, 3, 0, 0, 0.0) - Mjday(2025, 4, 3)) < TOL_);
    return 0;
}
int Mjday_02(){
    double x = Mjday(2025,4,12,1,1,1);
    _assert(fabs(Mjday(2025,4,12,1,1,1)-60777.0423726854)<TOL2_);
    return 0;
}
int R_x_01(){
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_x(alpha);

    _assert(fabs(sol(1,1) - 1 < TOL_ &&
    fabs(sol(1,2)) < TOL_ && fabs(sol(1,3)<TOL_)));
    _assert(sol(2,1)<TOL_ && fabs(sol(2,2)-0.54030230586814) < TOL_ && fabs(sol(2,3)
    -0.841470984807897) < TOL_);
    _assert(fabs(sol(1,1))- 1 < TOL_);

    return 0;

}
int R_y_01(){
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_y(alpha);

    _assert(fabs(sol(1,1) - 0.54030230586814) < TOL2_ && fabs(sol(1,3) + 0.841470984807897) < TOL2_);

    _assert(fabs(sol(2,1)-0<TOL2_) && fabs(sol(2,2)-1)<TOL2_ && fabs(sol(2,3)-0<TOL2_));

    _assert(fabs(sol(3,1)-0.841470984807897)<TOL2_ && fabs(sol(3,3)-0.54030230586814)<TOL2_);
    return 0;

}
int R_z_01(){
    double alpha = 1.0;
    Matrix sol(3,3);

    sol = R_z(alpha);

    _assert(fabs(sol(1,1) - 0.54030230586814) < TOL2_ && fabs(sol(1,2) - 0.841470984807897) < TOL2_);

    _assert(fabs(sol(2,1)+0.841470984807897<TOL2_) && fabs(sol(2,2)-0.54030230586814)<TOL2_ && fabs(sol(2,3)-0<TOL2_));

    _assert(fabs(sol(3,1)-0)<TOL2_ && fabs(sol(3,3)-1)<TOL2_);
    return 0;

}
int Legendre_01(){
    Matrix pnm(3,3), dpnm(3,3);
    Legendre(2,2,1,pnm,dpnm);
    pnm.print();
    dpnm.print();
    _assert(fabs(pnm(1,1)-1)<TOL_ && fabs(pnm(2,1)-1.4574704987823)<TOL2_ && fabs(pnm(3,1)-1.25691645573063)<TOL2_);
    _assert(fabs(pnm(1,2)-0)<TOL_ && fabs(pnm(2,2)-0.93583104521024)<TOL_ && fabs(pnm(3,2)-1.76084689542256)<TOL_);
    _assert(fabs(pnm(1,3)-0)<TOL_ && fabs(pnm(2,3)-0)<TOL_ && fabs(pnm(3,3)-0.56531339467086)<TOL_);

    _assert(fabs(dpnm(1,1)-0)<TOL_ && fabs(dpnm(2,1)-0.93583104521024)<TOL2_ && fabs(dpnm(3,1)-3.04987628722180)<TOL2_);
    _assert(fabs(dpnm(1,2)-0)<TOL_ && fabs(dpnm(2,2)+1.45747049878230)<TOL_ && fabs(dpnm(3,2)+1.61172976752398)<TOL_);
    _assert(fabs(dpnm(1,3)-0)<TOL_ && fabs(dpnm(2,3)-0)<TOL_ && fabs(dpnm(3,3)+1.76084689542256)<TOL_);

    return 0;
}
int SAT_Const(){
    _assert(fabs(consts.pi2 - 2 * M_PI) < TOL_);
    _assert(fabs(consts.Rad - M_PI / 180) < TOL_);
    _assert(fabs(consts.Deg - 180 / M_PI) < TOL_);
    _assert(fabs(consts.Arcs - 3600 * 180 / M_PI) < TOL_);

    _assert(fabs(consts.MJD_J2000 - 51544.5) < TOL_);
    _assert(fabs(consts.T_B1950 + 0.500002108) < TOL_);
    _assert(fabs(consts.c_light - 299792458.0) < TOL_);
    _assert(fabs(consts.AU - 149597870700.0) < TOL_);

    _assert(fabs(consts.R_Earth - 6378.1363e3) < TOL_);
    _assert(fabs(consts.f_Earth - 1 / 298.257223563) < TOL_);
    _assert(fabs(consts.R_Sun - 696000e3) < TOL_);
    _assert(fabs(consts.R_Moon - 1738e3) < TOL_);

    _assert(fabs(consts.omega_Earth - (15.04106717866910 / 3600 * (M_PI / 180))) < TOL_);

    _assert(fabs(consts.GM_Earth - 398600.435436e9) < TOL_);
    _assert(fabs(consts.GM_Sun - 132712440041.939400e9) < TOL_);
    _assert(fabs(consts.GM_Moon - (398600.435436e9 / 81.30056907419062)) < TOL_);
    _assert(fabs(consts.GM_Mercury - 22031.780000e9) < TOL_);
    _assert(fabs(consts.GM_Venus - 324858.592000e9) < TOL_);
    _assert(fabs(consts.GM_Mars - 42828.375214e9) < TOL_);
    _assert(fabs(consts.GM_Jupiter - 126712764.800000e9) < TOL_);
    _assert(fabs(consts.GM_Saturn - 37940585.200000e9) < TOL_);
    _assert(fabs(consts.GM_Uranus - 5794548.600000e9) < TOL_);
    _assert(fabs(consts.GM_Neptune - 6836527.100580e9) < TOL_);
    _assert(fabs(consts.GM_Pluto - 977.0000000000009e9) < TOL_);

    _assert(fabs(consts.P_Sol - 1367 / 299792458.0) < TOL_);

    return 0;
}
int AuxParam01() {
    _assert(fabs(auxParam.Mjd_UTC-4.974611635416653e+04)<TOL2_);
    _assert(fabs(auxParam.n-20)<TOL2_);
    _assert(fabs(auxParam.m-20)<TOL2_);
    _assert(fabs(auxParam.sun-1)<TOL2_);
    _assert(fabs(auxParam.moon-1)<TOL2_);
    _assert(fabs(auxParam.planets-1)<TOL2_);
    _assert(fabs(auxParam.Mjd_TT-4.974611706231468e+04)<TOL2_);
    return 0;
}
int NutAngles01(){
    double a,b;
    NutAngles(1123.567,a,b);
    _assert(fabs(a-7.99413796463657e-05)<TOL2_ && (fabs(b-1.42585698111339e-06)<TOL2_));

    return 0;
}
int Unit01(){
    double vec[] = {1,2,3};
    double outvec[3];
    unit(vec,outvec);
    _assert(fabs(outvec[0]-0.267261241912424)<TOL2_ && fabs(outvec[1]-0.534522483824849)<TOL2_
    && fabs(outvec[2]-0.801783725737273)<TOL2_);
    return 0;
}
int Unit02() {
    double vec[] = {1e-8,1e-8,1e-8};
    double outvec[3];
    unit(vec, outvec);
    _assert(outvec[0] == 0.0);
    _assert(outvec[1] == 0.0);
    _assert(outvec[2] == 0.0);
    return 0;
}
int AzElPa01(){
    double s[] ={1,2,3};
    double Az,El;
    double dAds[3], dEds[3];
    AzElPa(s,Az,El,dAds,dEds);
    _assert(fabs(Az-0.463647609000806)<TOL_ && fabs(El-0.930274014115472)<TOL_);
    _assert(fabs(dAds[0]-0.4)<TOL_ && fabs(dAds[1]+0.2)<TOL_ && fabs(dAds[2]-0)<TOL_);
    _assert(fabs(dEds[0]+0.095831484749991)<TOL_ && fabs(dEds[1]+0.191662969499982)<TOL_ && fabs(dEds[2]-0.159719141249985)<TOL_);
    return 0;
}
int Frac01(){
    _assert(fabs(Frac(5.2)-0.2)<TOL_);
    _assert(fabs(Frac(5)-0)<TOL_);
    _assert(fabs(Frac(2.1)-0.1)<TOL_);
    return 0;
}
int AccelPointMass01(){
    Matrix r(3,1), s(3,1);
    r(1,1) = 7078e3; r(2,1) = 0; r(3,1) = 0;
    s(1,1) = 384400e3; s(2,1) = 0; s(3,1) = 0;
    double GM = 4.9048695e12;

    Matrix result = AccelPointMass(r, s, GM);
    result.print();
    _assert(fabs(result(1,1)-1.25720640996716e-06)<TOL2_);
    _assert(result(2,1)==0);
    _assert(result(3,1)==0);

    return 0;
}
int IERS01(){
        double datos[] = {
                2025, 1, 1, 59780, 0.12, -0.15, -0.24, 0.0001, 0.05, -0.03, 0.0002, -0.0001, 37,
                2025, 1, 2, 59780, 0.13, -0.16, -0.23, 0.0002, 0.06, -0.02, 0.0003, -0.0002, 37,
                2025, 1, 3, 59782, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37,
                2025, 1, 3, 59780, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37,
                2025, 1, 3, 59780, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37,
                2025, 1, 3, 59780, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37,
                2025, 1, 3, 59780, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37,
                2025, 1, 3, 59780, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37,
                2025, 1, 3, 59780, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37,
                2025, 1, 3, 59780, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37,
                2025, 1, 3, 59780, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37,
                2025, 1, 3, 59780, 0.14, -0.17, -0.22, 0.0003, 0.07, -0.01, 0.0004, -0.0003, 37
        };
        Matrix m(13,13, datos,169);
        m = m.transpose();
        double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
        IERS(m, 59780.5, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, 'l');
        _assert(fabs(x_pole-0.0000006060)<TOL2_);
        _assert(fabs(y_pole+0.0000007515)<TOL2_);
        _assert(fabs(UT1_UTC+0.2350000000)<TOL2_);
        _assert(fabs(LOD-0.0001500000)<TOL2_);
        _assert(fabs(dpsi-0.0000002666)<TOL2_);
        _assert(fabs(deps+0.0000001212)<TOL2_);
        _assert(fabs(dx_pole-0.0000000012)<TOL2_);
        _assert(fabs(dy_pole+0.0000000007)<TOL2_);
        _assert(fabs(TAI_UTC-37.0)<TOL2_);
        return 0;
}
int MeanObliquity01(){
    _assert(fabs(MeanObliquity(1234567890)-339374.3819176)<TOL2_);
    return 0;
}

int TimeUpdate01() {
    double pValores[]={1,2,3,4};
    Matrix* P= new Matrix(2,2, pValores, 4);
    P->print();

    double phiValores[] = {5,6,7,8};
    Matrix* Phi= new Matrix(2,2,phiValores, 4);
    Phi->print();

    TimeUpdate(P,Phi);

    P->print();

    _assert(fabs((*P)(1,1)-319)<TOL_);
    _assert(fabs((*P)(1,2)-433)<TOL_);
    _assert(fabs((*P)(2,1)-431)<TOL_);
    _assert(fabs((*P)(2,2)-585)<TOL_);

    return 0;

}
int Mjday_TBD01() {
    double date = Mjday(2025,05,02,13,52,42);
    double x = Mjday_TDB(date);
    _assert(fabs(x-60797.57826390586706)<TOL2_);
    return 0;

}
int Identity01() {
    Matrix A(3,3);
    Matrix I = A.identity();

    for (int i=1;i<=3;i++) {
        for (int j=1; j<=3; j++) {
            if (i==j) {
                _assert(fabs(I(i,i)-1) < TOL_);
            }
            else {
                _assert(I(i,j)==0);
            }
        }
    }
    I.print();
    return 0;
}
int Inverse01() {
    double valores[] = {4, 7, 2, 6};
    Matrix A(2, 2, valores, 4);
    Matrix Ainv = A.inverse();

    _assert(fabs(Ainv(1, 1) - 0.6) < TOL_);
    _assert(fabs(Ainv(1, 2) + 0.7) < TOL_);
    _assert(fabs(Ainv(2, 1) + 0.2) < TOL_);
    _assert(fabs(Ainv(2, 2) - 0.4) < TOL_);

    Matrix I = A * Ainv;
    for (int i = 1; i <= 2; i++) {
        for (int j = 1; j <= 2; j++) {
            if (i == j) {
                _assert(fabs(I(i, j) - 1) < TOL_);
            } else {
                _assert(fabs(I(i, j)) < TOL_);
            }
        }
    }
    return 0;
}
int MeasUpdate01() {
    double valoresX[] = {1,1};
    Matrix x{2,1, valoresX, 2};
    cout << "Matrix x" << endl;
    x.print();
    double valoresZ[]{2,2};
    Matrix z{2,1, valoresZ, 2};
    cout << "Matrix z" << endl;
    z.print();
    double valoresg[] = {1,1};
    Matrix g{2,1,valoresg,2};
    cout << "Matrix g" << endl;
    g.print();
    double valoresS[] = {1,1};
    Matrix s{2,1, valoresS, 2};
    cout << "Matrix s" << endl;
    s.print();
    double valoresG[] = {1,0,0,1};
    Matrix G{2,2,valoresG, 4};
    cout << "Matrix G" << endl;
    G.print();
    double valoresP[] = {1,0,0,1};
    Matrix P{2,2,valoresP,4};
    cout << "Matrix P" << endl;
    P.print();
    Matrix K(2,2);
    cout << "Matrix K" << endl;
    K.print();
    MeasUpdate(x,z,g,s,G,P,2,K);

    _assert(fabs(K(1, 1) - 0.5) < TOL_);
    _assert(fabs(K(1, 2) - 0.0) < TOL_);
    _assert(fabs(K(2, 1) - 0.0) < TOL_);
    _assert(fabs(K(2, 2) - 0.5) < TOL_);


    // Verificación del vector x (2x1)
    _assert(fabs(x(1, 1) - 1.5) < TOL_);
    _assert(fabs(x(2, 1) - 1.5) < TOL_);

    // Verificación de la matriz P (2x2)
    _assert(fabs(P(1, 1) - 0.5) < TOL_);
    _assert(fabs(P(1, 2) - 0.0) < TOL_);
    _assert(fabs(P(2, 1) - 0.0) < TOL_);
    _assert(fabs(P(2, 2) - 0.5) < TOL_);

    return 0;
}

int Cheb3D01(){
    const int N = 3;
    const double Ta = 0.0;
    const double Tb = 1.0;
    const double t = 0.5;
    double Cx[] = {1.0, 0.5, -0.2};
    double Cy[] = {0.0, 1.0, 0.1};
    double Cz[] = {2.0, -0.3, 0.4};
    double ChebApp[3];
    Cheb3D(t,N,Ta,Tb,Cx,Cy,Cz,ChebApp);
    _assert(fabs(ChebApp[0]-1.2)<TOL_);
    _assert(fabs(ChebApp[1]+0.1)<TOL_);
    _assert(fabs(ChebApp[2]-1.6)<TOL_);
    return 0;
}
int Sign_01() {

    _assert(sign_(5,3) == 5);
    _assert(sign_(5,-3) == -5);
    _assert(sign_(-5,-3)== -5);
    _assert(sign_(4.2,0) == 4.2);
    _assert(sign_(-4.2,0) == 4.2);

    return 0;
}
int NutMatrix01(){
    Matrix x = NutMatrix(12345678);
    _assert(fabs(x(1,1)-0.999999996289796)<TOL2_ );
    _assert(fabs(x(1,2)+6.76497562849557e-05)<TOL2_);
    _assert(fabs(x(1,3)+5.33284056190134e-05)<TOL2_);
    _assert(fabs(x(2,1)-6.76483187988301e-05)<TOL2_);
    _assert(fabs(x(2,2)-0.999999997348521)<TOL2_);
    _assert(fabs(x(2,3)+2.6956698330638e-05)<TOL2_);
    _assert(fabs(x(3,1)-5.33302290916865e-05)<TOL2_);
    _assert(fabs(x(3,2)-2.69530906537585e-05)<TOL2_);
    _assert(fabs(x(3,3)-0.999999998214709)<TOL2_);

    return 0;
}
int PoleMatrix01(){
    Matrix x = PoleMatrix(12.4,48.2);
    _assert(fabs(x(1,1)-0.986192302278864)<TOL2_);
    _assert(fabs(x(1,2)-0.145751176665682)<TOL2_);
    _assert(fabs(x(1,3)-0.078621481965705)<TOL2_);
    _assert(x(2,1)==0);
    _assert(fabs(x(2,2)+0.474755432662659)<TOL2_);
    _assert(fabs(x(2,3)-0.880117764368662)<TOL2_);
    _assert(fabs(x(3,1)-0.165604175448309)<TOL2_);
    _assert(fabs(x(3,2)+0.867965364319257)<TOL2_);
    _assert(fabs(x(3,3)+0.468200153156986)<TOL2_);

    return 0;
}

int PrecMatrix01(){
    Matrix x = PrecMatrix(12345678,87654321);
    _assert(fabs(x(1,1)-0.817653102598568)<TOL2_);
    _assert(fabs(x(1,2)-0.462267643134377)<TOL2_);
    _assert(fabs(x(1,3)+0.343150156523239)<TOL2_);
    _assert(fabs(x(2,1)+0.567591358508052)<TOL2_);
    _assert(fabs(x(2,2)-0.547521381081051)<TOL2_);
    _assert(fabs(x(2,3)+0.614866153732731)<TOL2_);
    _assert(fabs(x(3,1)+0.0963506801113464)<TOL2_);
    _assert(fabs(x(3,2)-0.697516281795691)<TOL2_);
    _assert(fabs(x(3,3)-0.710061675540931)<TOL2_);

    return 0;
}
int timediff01(){
    double a,b,c,d,e;
    timediff(24,24,a,b,c,d,e);
    _assert(a==0);
    _assert(fabs(b+5)<TOL2_);
    _assert(fabs(c-19)<TOL2_);
    _assert(fabs(d-56.184)<TOL2_);
    _assert(fabs(e-5)<TOL2_);

    return 0;
}
int EccAnom01(){
    double x = EccAnom(135,0.9);
    _assert(fabs(x-3.09501409518025)<TOL_);
    return 0;
}
int elements01(){
    double p,a,e,i,Omega,omega,M;
    double valores[] = {7000e3, 0,0,0, 7.5e3,0};
    elements(valores,p,a,e,i,Omega,omega,M);
    _assert(fabs(p-6914819.34028782)<TOL2_);

    _assert(fabs(a-6915843.4136513)<TOL3_);

    _assert(fabs(e-0.0121686656731681)<TOL2_);
    _assert(i==0 && Omega==0);
    _assert(fabs(omega-3.14159265358979)<TOL2_);
    _assert(fabs(M-3.14159265358979)<TOL2_);
    return 0;
    /*
     * Valores para el test:
     * p=6914819.34028782
     * a=6915843.4136513
     * e=0.0121686656731681
     * i=0
     * Omega=0
     * omega=3.14159265358979
     * M=3.14159265358979
     */
}
int EqnEquinox01(){
    double res = EqnEquinox(1234567890);
    _assert(fabs(res+0.00169270155393772)<TOL2_);
    return 0;
}
int gmst01(){
    _assert(fabs(gmst(1234567890)-2.51251437041138)<TOL2_);
    return 0;
}
int gast01(){
    _assert(fabs(gast(1234567890)-2.51082166885744));
    return 0;
}
int readFromFile01(){
    Matrix mat = Matrix::LoadFromFile("eop19620101.txt");

    Matrix matTras = Matrix(mat.getCol(), mat.getFil());
    matTras = mat.transpose();
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    IERS(matTras, 37668, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, 'l');

    _assert(fabs(x_pole+1.06654161707287e-07)<TOL_);
    _assert(fabs(y_pole-1.04865684037674e-06)<TOL_);
    _assert(fabs(UT1_UTC-0.0311435)<TOL_);
    _assert(fabs(LOD-0.001496)<TOL_);
    _assert(fabs(dpsi-3.09762005271316e-07)<TOL_);
    _assert(fabs(deps-3.18716513961409e-08)<TOL_);
    _assert(fabs(TAI_UTC-2)<TOL_);
    _assert(fabs(dx_pole-0)<TOL_);
    _assert(fabs(dy_pole-0)<TOL_);

    return 0;
}
int GHAMatrix01() {
    Matrix mat = GHAMatrix(1234567890);
    _assert(fabs(mat(1,1)+0.807573046590817)<TOL2_);
    _assert(fabs(mat(1,2)-0.589767559653823)<TOL2_);
    _assert(fabs(mat(1,3)-0)<TOL2_);
    _assert(fabs(mat(2,1)+0.589767559653823)<TOL2_);
    _assert(fabs(mat(2,2)+0.807573046590817)<TOL2_);
    _assert(fabs(mat(2,3)-0)<TOL2_);
    _assert(fabs(mat(3,1)-0)<TOL2_);
    _assert(fabs(mat(3,2)-0)<TOL2_);
    _assert(fabs(mat(3,3)-1)<TOL2_);
    return 0;
}
int JPL_EphDE43001() {
    Matrix r_Mercury(3, 1);
    Matrix r_Venus(3, 1);
    Matrix r_Earth(3, 1);
    Matrix r_Mars(3, 1);
    Matrix r_Jupiter(3, 1);
    Matrix r_Saturn(3, 1);
    Matrix r_Uranus(3, 1);
    Matrix r_Neptune(3, 1);
    Matrix r_Pluto(3, 1);
    Matrix r_Moon(3, 1);
    Matrix r_Sun(3, 1);

    JPL_Eph_DE430(58849.5, r_Earth, r_Mars, r_Mercury, r_Venus, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

    //Mercury
    _assert(fabs(r_Mercury(1,1)-18364570296.7767)<TOL4_);
    _assert(fabs(r_Mercury(2,1)+194240890308.454)<TOL4_);
    _assert(fabs(r_Mercury(3,1)+89580185100.4172)<TOL4_);

    //Venus
    _assert(fabs(r_Venus(1,1) - 134236480603.468) < TOL4_);
    _assert(fabs(r_Venus(2,1) - -121788430081.357) < TOL4_);
    _assert(fabs(r_Venus(3,1) - -59451902859.1386) < TOL4_);

    // Earth
    _assert(fabs(r_Earth(1,1) - -26742322957.984) < TOL4_);
    _assert(fabs(r_Earth(2,1) - 133827327107.372) < TOL4_);
    _assert(fabs(r_Earth(3,1) - 58018326578.9664) < TOL4_);

    // Mars
    _assert(fabs(r_Mars(1,1) - -170687646029.87) < TOL4_);
    _assert(fabs(r_Mars(2,1) - -255905377504.909) < TOL4_);
    _assert(fabs(r_Mars(3,1) - -108721468283.956) < TOL4_);

    // Jupiter
    _assert(fabs(r_Jupiter(1,1) - 105439249049.906) < TOL4_);
    _assert(fabs(r_Jupiter(2,1) - -847168695382.084) < TOL4_);
    _assert(fabs(r_Jupiter(3,1) - -365696885445.354) < TOL4_);

    // Saturn
    _assert(fabs(r_Saturn(1,1) - 594596636534.412) < TOL4_);
    _assert(fabs(r_Saturn(2,1) - -1408094254368.45) < TOL4_);
    _assert(fabs(r_Saturn(3,1) - -608810768966.399) < TOL4_);

    // Uranus
    _assert(fabs(r_Uranus(1,1) - 2453302678658.18) < TOL4_);
    _assert(fabs(r_Uranus(2,1) - 1439200900737.1) < TOL4_);
    _assert(fabs(r_Uranus(3,1) - 596605466768.131) < TOL4_);

    // Neptune
    _assert(fabs(r_Neptune(1,1) - 4400884479114.87) < TOL4_);
    _assert(fabs(r_Neptune(2,1) - -974205297187.0) < TOL4_);
    _assert(fabs(r_Neptune(3,1) - -510890407765.413) < TOL4_);

    // Pluto
    _assert(fabs(r_Pluto(1,1) - 1967686040254.54) < TOL4_);
    _assert(fabs(r_Pluto(2,1) - -4414806220850.92) < TOL4_);
    _assert(fabs(r_Pluto(3,1) - -1978780908608.13) < TOL4_);

    // Moon
    _assert(fabs(r_Moon(1,1) - 398673022.292818) < TOL4_);
    _assert(fabs(r_Moon(2,1) - -38480993.8201334) < TOL4_);
    _assert(fabs(r_Moon(3,1) - -55662807.01252) < TOL4_);

    // Sun
    _assert(fabs(r_Sun(1,1) - 26173432883.2097) < TOL4_);
    _assert(fabs(r_Sun(2,1) - -132807686289.682) < TOL4_);
    _assert(fabs(r_Sun(3,1) - -57572484281.2736) < TOL4_);
    return 0;
}
int LTC01() {
    Matrix M = LTC(0.5,0.75);

    _assert(fabs(M(1,1)+0.479425538604203)<TOL2_);
    _assert(fabs(M(1,2)-0.877582561890373)<TOL2_);
    _assert(fabs(M(1,3)-0.0)<TOL2_);
    _assert(fabs(M(2,1)+0.598194289305055)<TOL2_);
    _assert(fabs(M(2,2)+0.326795029657688)<TOL2_);
    _assert(fabs(M(2,3)-0.731688868873821)<TOL2_);
    _assert(fabs(M(3,1)-0.642117392052957)<TOL2_);
    _assert(fabs(M(3,2)-0.350790330050532)<TOL2_);
    _assert(fabs(M(3,3)-0.681638760023334)<TOL2_);

    return 0;

}
int Cnm01() {
    _assert(fabs(Cnm(181,181)-1.52619566516600e-09)<TOL2_);
    _assert(fabs(Cnm(179,81)+1.38673296550100e-09)<TOL2_);
    return 0;
}
int Snm01(){
    _assert(fabs(Snm(181,181)+1.99213234732500e-09)<TOL2_);
    _assert(fabs(Snm(179,81)+3.43860639511800e-10)<TOL2_);
    return 0;
}
int obs01() {
    _assert(fabs(obs(46,4)-2653472)<TOL2_);
    _assert(fabs(obs(42,1)-49746.1158796297)<TOL2_);
    return 0;
}
int accelHarmonic01() {
    double valores[]= {6378136,0,0};
    Matrix r(3,1,valores,3);
    int n_max = 4, m_max = n_max;

    Matrix E1 (3,3);
    Matrix E = E1.identity();
    Matrix res = AccelHarmonic(r,E,n_max,m_max);
    _assert(fabs(res(1,1)+9.81425003387253)<TOL2_);
    _assert(fabs(res(2,1)-3.349605551004308e-05)<TOL4_);
    _assert(fabs(res(3,1)-0.000106491638754056)<TOL2_);

    return 0;
}

int GaccelHarmonic01() {
    double valores[]= {7000e3,0,0};
    Matrix r(3,1,valores,3);
    Matrix U1 (3,3);
    Matrix U = U1.identity();

    Matrix G = G_AccelHarmonic(r,U,4,4);

    _assert(fabs(G(1,1)-2.33048175246608e-06)<TOL2_);
    _assert(fabs(G(2,1)+1.88080662155699e-11)<TOL2_);
    _assert(fabs(G(3,1)+5.18909359925601e-11)<TOL2_);

    _assert(fabs(G(1,2)+1.88040301950557e-11)<TOL2_);
    _assert(fabs(G(2,2)+1.1636930401132e-06)<TOL2_);
    _assert(fabs(G(1,2)+5.90145164317381e-12)<TOL2_);

    _assert(fabs(G(1,3)+5.18900588865551e-11)<TOL2_);
    _assert(fabs(G(2,3)+5.90145164995007e-12)<TOL2_);
    _assert(fabs(G(3,3)+1.16678871423606e-06)<TOL2_);

    return 0;
}
int VarEqn01() {
    Matrix yPhi(42,1);
    yPhi(1,1) = 1e6;
    yPhi(2,1) = 2e6;
    yPhi(3,1) = 3e6;
    yPhi(4,1) = 1e3;
    yPhi(5,1) = 2e3;
    yPhi(6,1) = 3e3;
    for (int i = 7; i <= 42; i++) {
        yPhi(i,1) = 0.0;
    }

    Matrix result = VarEqn(4,yPhi);

    _assert(fabs(result(1,1)-1000)<TOL2_);
    _assert(fabs(result(2,1)-2000)<TOL2_);
    _assert(fabs(result(3,1)-3000)<TOL2_);
    _assert(fabs(result(4,1)+8.55483473704505)<TOL4_);
    _assert(fabs(result(5,1)+14.7125328326)<TOL2_);
    _assert(fabs(result(6,1)+26.64751645926)<TOL2_);

    return 0;
}
int Accel01() {
    cout << "Accel" << endl;
    double valores[] = {4.6,3.8,0,9.0,7.5,2.8};
    Matrix A(6,1, valores, 6);

    Matrix B = Accel(0,A);
    B.print();
    
    _assert(fabs(B(1,1)-9)<TOL2_);
    _assert(fabs(B(2,1)-7.5)<TOL2_);
    _assert(fabs(B(3,1)-2.8)<TOL2_);


    return 0;
}

int all_tests()
{
    _verify(testConstructorVacio);
    _verify(testConstructorVector);
    _verify(testSuma);
    _verify(testResta);
    _verify(testProducto);
    _verify(testTraspuesta);
    _verify(Identity01);
    _verify(Inverse01);
    _verify(Mjday_01);
    _verify(Mjday_02);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(Legendre_01);
    _verify(SAT_Const);
    _verify(NutAngles01);
    _verify(Unit01);
    _verify(Unit02);
    _verify(AzElPa01);
    _verify(Frac01);
    _verify(IERS01);
    _verify(MeanObliquity01);
    _verify(TimeUpdate01);
    _verify(Mjday_TBD01);
    _verify(MeasUpdate01);
    _verify(Cheb3D01);
    _verify(Sign_01);
    _verify(NutMatrix01);
    _verify(PoleMatrix01);
    _verify(PrecMatrix01);
    _verify(timediff01);
    _verify(EccAnom01);
    _verify(elements01);
    _verify(EqnEquinox01);
    _verify(gmst01);
    _verify(gast01);
    _verify(readFromFile01);
    _verify(AuxParam01);
    //_verify(JPL_EphDE43001);
    _verify(GHAMatrix01);
    _verify(LTC01);
    _verify(Cnm01);
    _verify(Snm01);
    _verify(obs01);
    _verify(accelHarmonic01);
    _verify(GaccelHarmonic01);
    _verify(VarEqn01);
    //_verify(AccelPointMass01);
    _verify(Accel01);

    return 0;
}


int main()
{
    GGM03S();
    GEOS3();
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");
    printf("Tests run: %d\n", tests_run);

    return result != 0;
}


