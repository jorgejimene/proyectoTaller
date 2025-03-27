/*
 * Archivo de tests
 *
 */

#include <stdio.h>
#include "./include/Matrix.h"

int tests_run = 0;

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


int all_tests()
{
    _verify(testConstructorVacio);
    _verify(testConstructorVector);
    _verify(testSuma);
    _verify(testResta);
    _verify(testProducto);
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


