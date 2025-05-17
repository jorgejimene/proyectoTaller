#include <Matrix.h>

#include "Cheb3D.h"
/*
%--------------------------------------------------------------------------
%
% JPL_Eph_DE430: Computes the sun, moon, and nine major planets' equatorial
%                position using JPL Ephemerides
%
% Inputs:
%   Mjd_TDB         Modified julian date of TDB
%
% Output:
%   r_Earth(solar system barycenter (SSB)),r_Mars,r_Mercury,r_Venus,
%   r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,
%   r_Sun(geocentric equatorial position ([m]) referred to the
%   International Celestial Reference Frame (ICRF))
%
% Notes: Light-time is already taken into account
%
% Last modified:   2018/01/11   M. Mahooti
%
%--------------------------------------------------------------------------
*/
void JPL_Eph_DE430(double Mjd_TDB,double r_Mercury[3],double r_Venus[3],double r_Earth[3],double r_Mars[3],double r_Jupiter[3],double r_Saturn[3],double r_Uranus[3], double r_Neptune[3],double r_Pluto[3],double r_Moon[3],double r_Sun[3]) {
    Matrix PC = Matrix::LoadFromFile("PC.txt");
    int rowCount = PC.getFil();
    int colCount = PC.getCol();
    int i=-1;
    double JD = Mjd_TDB + 2400000.5;

    for (int il = 0; il < rowCount; ++il) {
        // Accedemos a la primera y segunda columna (indexada desde 0 en C++)
        double col1 = PC(il + 1, 1); // Primera columna (ajustar a índice 1 en Matrix)
        double col2 = PC(il + 1, 2); // Segunda columna (ajustar a índice 1 en Matrix)

        // Verificar si JD está entre col1 y col2
        if (col1 <= JD && JD <= col2) {
            i=il+1;// Encontramos el índice
            break;
        }
    }
    double PCtemp[colCount];
    for (int j=0;j<=colCount; j++) {
        PCtemp[j] = PC(i,j+1);
    }
    double t1 = PCtemp[0] - 2400000.5;
    double dt = Mjd_TDB - t1;


    double EMRAT = 81.30056907419062;
    double EMRAT1 = 1/(1+EMRAT);

    for (int k=0; k<3; k++) {
        r_Earth[k] -= EMRAT1 * r_Moon[k];
    }
    for (int k=0; k<3; k++) {
        r_Mercury[k] -= r_Earth[k];
        r_Venus[k] -= r_Earth[k];
        r_Mars[k] -= r_Earth[k];
        r_Jupiter[k] -= r_Earth[k];
        r_Saturn[k] -= r_Earth[k];
        r_Uranus[k] -= r_Earth[k];
        r_Neptune[k] -= r_Earth[k];
        r_Pluto[k] -= r_Earth[k];
        r_Sun[k] -= r_Earth[k];
    }
}
