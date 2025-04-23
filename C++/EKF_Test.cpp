/*
 * Archivo de tests
 *
 */

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "./INCLUDE/Matrix.h"
#include "./INCLUDE/Mjday.h"
#include "./INCLUDE/R_x.h"
#include "./INCLUDE/R_y.h"
#include "./INCLUDE/R_z.h"
#include "./INCLUDE/Legendre.h"
#include "./INCLUDE/SAT_Const.h"
#include "./INCLUDE/NutAngles.h"
#include "./INCLUDE/Unit.h"

int tests_run = 0;

#define TOL_ 10e-14
#define TOL2_ 10e-10 //para numeros mas pequeños
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
int all_tests()
{
    _verify(testConstructorVacio);
    _verify(testConstructorVector);
    _verify(testSuma);
    _verify(testResta);
    _verify(testProducto);
    _verify(Mjday_01);
    _verify(Mjday_02);
    _verify(R_x_01);
    _verify(R_y_01);
    _verify(R_z_01);
    _verify(Legendre_01);
    _verify(SAT_Const);
    _verify(NutAngles01);
    _verify(Unit01);

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}


