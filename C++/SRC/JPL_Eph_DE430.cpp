#include "../INCLUDE/JPL_Eph_DE430.h"



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
void JPL_Eph_DE430(double Mjd_TDB, Matrix& r_Earth, Matrix& r_Mars, Matrix& r_Mercury, Matrix& r_Venus,
                   Matrix& r_Jupiter, Matrix& r_Saturn, Matrix& r_Uranus, Matrix& r_Neptune,
                   Matrix& r_Pluto, Matrix& r_Moon, Matrix& r_Sun) {
    Matrix PC = Matrix::LoadFromFile("PC.txt");
    double JD = Mjd_TDB + 2400000.5;

    int rowCount = PC.getFil();
    int colCount = PC.getCol();
    int i=0;
    double EMRAT = 81.30056907419062;
    double EMRAT1 = 1/(1+EMRAT);

    for (int il = 1; il <= rowCount; ++il) {
        if (PC(il,1) <= JD && JD <= PC(il,2)) {
            i=il;
            break;
        }
    }
    Matrix PCtemp(1,colCount);
    for (int j=1;j<=colCount; j++) {
        PCtemp(1,j) = PC(i,j);
    }
    double t1 = PCtemp(1,1) - 2400000.5;
    double dt = Mjd_TDB - t1;

    //Cx_Earth = PCtemp(temp(1):temp(2)-1);
    //Cy_Earth = PCtemp(temp(2):temp(3)-1);
    //Cz_Earth = PCtemp(temp(3):temp(4)-1);

    //Earth
    Matrix Cx_Earth(26,1), Cy_Earth(26,1), Cz_Earth(26,1);
    for (int c = 0; c < 13; c++) {
        Cx_Earth(c+1,1) =PCtemp(1,231+c);
        Cy_Earth(c+1,1) = PCtemp(1,244+c);
        Cz_Earth(c+1,1) = PCtemp(1,257+c);
    }
    //temp = temp+39;
    //Cx = PCtemp(temp(1):temp(2) - 1);
    //Cy = PCtemp(temp(2):temp(3) - 1);
    //Cz = PCtemp(temp(3):temp(4) - 1);
    //Cx_Earth = [Cx_Earth,Cx];
    //Cy_Earth = [Cy_Earth,Cy];
    //Cz_Earth = [Cz_Earth,Cz];
    for (int c = 0; c < 13; c++) {
        Cx_Earth(c+14,1) = PCtemp(1,270+c);
        Cy_Earth(c+14,1) = PCtemp(1,283+c);
        Cz_Earth(c+14,1) = PCtemp(1,296+c);
    }
    int j;
    double Mjd0;
    if (0<= dt && dt <= 16) {
        j=0;
        Mjd0 = t1;
    }
    else if (16<dt && dt <=32) {
        j=1;
        Mjd0 = t1+16*j;
    }

    double Cx_Earth2[13],Cy_Earth2[13],Cz_Earth2[13];
    for (int c = 0; c < 13; c++) {
        Cx_Earth2[c] = Cx_Earth(j*13+c+1,1);
        Cy_Earth2[c]= Cy_Earth(j*13+c+1,1);
        Cz_Earth2[c]= Cz_Earth(j*13+c+1,1);
    }


    double res[3];
    Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16,Cx_Earth2, Cy_Earth2,Cz_Earth2, res);
    for (int c = 0; c < 3; c++) {
        r_Earth(c+1,1) = res[c]*1e3;
    }

    //Moon
    Matrix Cx_Moon(104,1),Cy_Moon(104,1),Cz_Moon(104,1);
    for (int c = 0; c < 13; c++) {
        Cx_Moon(c+1,1) = PCtemp(1,441+c);
        Cy_Moon(c+1,1) = PCtemp(1,454+c);
        Cz_Moon(c+1,1) = PCtemp(1,467+c);
    }

    for (int m=1; m<=7;++m) {
        for (int c=0;c<13;++c) {
            Cx_Moon(13 * m + c +1, 1)= PCtemp(1,441+39*m+c);
            Cy_Moon(13*m + c +1, 1)= PCtemp(1,454+39*m+c);
            Cz_Moon(13*m + c +1, 1)= PCtemp(1,467+39*m+c);
        }
    }
    if (0<=dt && dt<=4) {
        j=0;
        Mjd0 = t1;
    }
    else if(4<dt && dt<=8) {
        j=1;
        Mjd0 = t1+4*j;
    }
    else if(8<dt && dt<=12) {
        j=2;
        Mjd0 = t1+4*j;
    }
    else if(12<dt && dt<=16) {
        j=3;
        Mjd0 = t1+4*j;
    }
    else if(16<dt && dt<=20) {
        j=4;
        Mjd0 = t1+4*j;
    }
    else if(20<dt && dt<=24) {
        j=5;
        Mjd0 = t1+4*j;
    }
    else if(24<dt && dt<=28) {
        j=6;
        Mjd0 = t1+4*j;
    }
    else if(28<dt && dt<=32) {
        j=7;
        Mjd0 = t1+4*j;
    }

    double Cx_Moon2[13],Cy_Moon2[13],Cz_Moon2[13];
    for (int c = 0; c < 13; c++) {
        Cx_Moon2[c] = Cx_Moon(j*13+c+1,1);
        Cy_Moon2[c] = Cy_Moon(j*13+c+1,1);
        Cz_Moon2[c] = Cz_Moon(j*13+c+1,1);
    }
    Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4,Cx_Moon2, Cy_Moon2,Cz_Moon2, res);
    for (int c = 0; c < 3; c++) {
        r_Moon(c+1,1) = res[c]*1e3;
    }



    //Sun
    Matrix Cx_Sun(22,1), Cy_Sun(22,1),Cz_Sun(22,1);
    for (int c=0; c<11;++c) {
        Cx_Sun(c+1,1) = PCtemp(1,753+c);
        Cy_Sun(c+1,1) = PCtemp(1,764+c);
        Cz_Sun(c+1,1) = PCtemp(1,775+c);
    }
    for (int c=0; c<11;++c) {
        Cx_Sun(c+12,1) = PCtemp(1,786+c);
        Cy_Sun(c+12,1) = PCtemp(1,797+c);
        Cz_Sun(c+12,1) = PCtemp(1,808+c);
    }
    if (0<=dt && dt<=16) {
        j=0;
        Mjd0 = t1;
    }
    else if(16<dt && dt<=32) {
        j=1;
        Mjd0 = t1+16*j;
    }
    double Cx_Sun2[11],Cy_Sun2[11],Cz_Sun2[11];
    for (int c = 0; c < 11; c++) {
        Cx_Sun2[c] = Cx_Sun(j*11+c+1,1);
        Cy_Sun2[c] = Cy_Sun(j*11+c+1,1);
        Cz_Sun2[c] = Cz_Sun(j*11+c+1,1);
    }
    Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16,Cx_Sun2, Cy_Sun2,Cz_Sun2, res);
    for (int c = 0; c < 3; c++) {
        r_Sun(c+1,1) = res[c]*1e3;
    }

    //Mercury
    Matrix Cx_Mercury(56,1),Cy_Mercury(56,1),Cz_Mercury(56,1);
    for (int c = 0; c < 14; ++c) {
        Cx_Mercury(c+1,1) = PCtemp(1,3+c);
        Cy_Mercury(c+1,1) = PCtemp(1,17+c);
        Cz_Mercury(c+1,1) = PCtemp(1,31+c);
    }
    for (int m = 1; m<=3; ++m) {
        for (int c=0; c<14; c++) {
            Cx_Mercury(14*m+c+1,1) = PCtemp(1,3+42*m+c);
            Cy_Mercury(14*m+c+1,1) = PCtemp(1,17+42*m+c);
            Cz_Mercury(14*m+c+1,1) = PCtemp(1,31+42*m+c);
        }
    }
    if (0<=dt && dt<=8) {
        j=0;
        Mjd0 = t1;
    }
    else if(8<dt && dt<=16) {
        j=1;
        Mjd0 = t1+8*j;
    }
    else if (16<dt && dt<=24) {
        j=2;
        Mjd0 = t1+8*j;
    }
    else if(24<dt && dt<=32) {
        j=3;
        Mjd0 = t1+8*j;
    }
    double Cx_Mercury2[14],Cy_Mercury2[14],Cz_Mercury2[14];
    for (int c = 0; c < 14; ++c) {
        Cx_Mercury2[c] = Cx_Mercury(j*14+c+1,1);
        Cy_Mercury2[c] = Cy_Mercury(j*14+c+1,1);
        Cz_Mercury2[c] = Cz_Mercury(j*14+c+1,1);
    }

    Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8,Cx_Mercury2, Cy_Mercury2,Cz_Mercury2, res);
    for (int c = 0; c < 3; c++) {
        r_Mercury(c+1,1) = res[c]*1e3;
    }


    //Venus
    Matrix Cx_Venus(20,1),Cy_Venus(20,1),Cz_Venus(20,1);
    for (int c=0;c<10;++c) {
        Cx_Venus(c+1,1) = PCtemp(1,171+c);
        Cy_Venus(c+1,1) = PCtemp(1,181+c);
        Cz_Venus(c+1,1) = PCtemp(1,191+c);
    }
    for (int c=0; c<10; ++c) {
        Cx_Venus(10+c+1,1) = PCtemp(1,201+c);
        Cy_Venus(10+c+1,1) = PCtemp(1,211+c);
        Cz_Venus(10+c+1,1) = PCtemp(1,221+c);
    }
    if (0<=dt && dt<=16) {
        j=0;
        Mjd0 = t1;
    }
    else if(16<dt && dt<=32) {
        j=1;
        Mjd0 = t1+16*j;
    }
    double Cx_Venus2[10],Cy_Venus2[10],Cz_Venus2[10];
    for (int c = 0; c < 10; ++c) {
        Cx_Venus2[c] = Cx_Venus(j*10+c+1,1);
        Cy_Venus2[c] = Cy_Venus(j*10+c+1,1);
        Cz_Venus2[c] = Cz_Venus(j*10+c+1,1);
    }
    Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16,Cx_Venus2, Cy_Venus2,Cz_Venus2, res);
    for (int c = 0; c < 3; c++) {
        r_Venus(c+1,1) = res[c]*1e3;
    }


    //Mars

    double Cx_Mars2[11],Cy_Mars2[11],Cz_Mars2[11];
    for (int c = 0; c < 11; ++c) {
        Cx_Mars2[c]= PCtemp(1,309+c);
        Cy_Mars2[c] = PCtemp(1,320+c);
        Cz_Mars2[c] = PCtemp(1,331+c);
    }
    j=0;
    Mjd0 = t1;
    Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32,Cx_Mars2, Cy_Mars2,Cz_Mars2, res);
    for (int c = 0; c < 3; c++) {
        r_Mars(c+1,1) = res[c]*1e3;
    }


    //Jupiter
    double Cx_Jupiter[8],Cy_Jupiter[8],Cz_Jupiter[8];
    for (int c=0; c < 8; ++c) {
        Cx_Jupiter[c]= PCtemp(1,342+c);
        Cy_Jupiter[c] = PCtemp(1,350+c);
        Cz_Jupiter[c] = PCtemp(1,358+c);
    }
    j=0;
    Mjd0 = t1;
    Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32,Cx_Jupiter, Cy_Jupiter,Cz_Jupiter, res);
    for (int c = 0; c < 3; c++) {
        r_Jupiter(c+1,1) = res[c]*1e3;
    }



    //Saturn

    double Cx_Saturn[7],Cy_Saturn[7],Cz_Saturn[7];
    for (int c = 0; c < 7; ++c) {
        Cx_Saturn[c]= PCtemp(1,366+c);
        Cy_Saturn[c] = PCtemp(1,373+c);
        Cz_Saturn[c] = PCtemp(1,380+c);
    }
    j=0;
    Mjd0 = t1;
    Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32,Cx_Saturn, Cy_Saturn,Cz_Saturn, res);
    for (int c = 0; c < 3; c++) {
        r_Saturn(c+1,1) = res[c]*1e3;
    }


    //Uranus

    double Cx_Uranus[6],Cy_Uranus[6],Cz_Uranus[6];
    for (int c = 0; c < 6; ++c) {
        Cx_Uranus[c]= PCtemp(1,387+c);
        Cy_Uranus[c] = PCtemp(1,393+c);
        Cz_Uranus[c] = PCtemp(1,399+c);
    }
    j=0;
    Mjd0 = t1;
    Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32,Cx_Uranus, Cy_Uranus,Cz_Uranus, res);
    for (int c = 0; c < 3; c++) {
        r_Uranus(c+1,1) = res[c]*1e3;
    }


    //Neptune
    double Cx_Neptune[6],Cy_Neptune[6],Cz_Neptune[6];
    for (int c = 0; c < 6; ++c) {
        Cx_Neptune[c]= PCtemp(1,405+c);
        Cy_Neptune[c] = PCtemp(1,411+c);
        Cz_Neptune[c] = PCtemp(1,417+c);
    }
    j=0;
    Mjd0 = t1;
    Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32,Cx_Neptune, Cy_Neptune,Cz_Neptune, res);
    for (int c = 0; c < 3; c++) {
        r_Neptune(c+1,1) = res[c]*1e3;
    }




    //Pluto
    double Cx_Pluto[6],Cy_Pluto[6],Cz_Pluto[6];
    for (int c = 0; c < 6; ++c) {
        Cx_Pluto[c]= PCtemp(1,423+c);
        Cy_Pluto[c] = PCtemp(1,429+c);
        Cz_Pluto[c] = PCtemp(1,435+c);
    }
    j=0;
    Mjd0 = t1;
    Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32,Cx_Pluto, Cy_Pluto,Cz_Pluto, res);
    for (int c = 0; c < 3; c++) {
        r_Pluto(c+1,1) = res[c]*1e3;
    }

    //Nutations
    Matrix Cx_Nutations(40,1), Cy_Nutations(40,1);
    for (int c=0; c < 10; ++c) {
        Cx_Nutations(c+1,1) = PCtemp(1,819+c);
        Cy_Nutations(c+1,1) = PCtemp(1,829+c);
    }
    for (int m = 1; m<= 3; ++m) {
        for (int c = 0; c < 10; ++c) {
            Cx_Nutations(10*m+c+1,1) = PCtemp(1,819+20*m+c);
            Cy_Nutations(10*m+c+1,1) = PCtemp(1,829+20*m+c);
        }
    }
    if (0<=dt && dt<=8) {
        j=0;
        Mjd0 = t1;
    }
    else if(8<dt && dt<=16) {
        j=1;
        Mjd0 = t1+8*j;
    }
    else if (16<dt && dt<=24) {
        j=2;
        Mjd0 = t1+8*j;
    }
    else if(24<dt && dt<=32) {
        j=3;
        Mjd0 = t1+8*j;
    }
    double Cx_Nutations2[10],Cy_Nutations2[10];
    double Cz_Nutations2[]={0,0,0,0,0,0,0,0,0,0};
    for (int c = 0; c < 10; ++c) {
        Cx_Nutations2[c] = Cx_Nutations(j*10+c+1,1);
        Cy_Nutations2[c] = Cy_Nutations(j*10+c+1,1);
    }

    Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8,Cx_Nutations2, Cy_Nutations2,Cz_Nutations2, res);
    Matrix Nutations(3,1);
    for (int c = 0; c < 3; c++) {
        Nutations(c+1,1) = res[c];
    }


    //Librations
    Matrix Cx_Librations(40,1), Cy_Librations(40,1), Cz_Librations(40,1);
    for (int c = 0; c < 10; ++c) {
        Cx_Librations(c+1,1) = PCtemp(1,899+c);
        Cy_Librations(c+1,1) = PCtemp(1,909+c);
        Cz_Librations(c+1,1) = PCtemp(1,919+c);
    }
    for (int m=1; m<= 3; ++m) {
        for (int c = 0; c < 10; ++c) {
            Cx_Librations(10*m+c+1,1) = PCtemp(1,899+30*m+c);
            Cy_Librations(10*m+c+1,1) = PCtemp(1,909+30*m+c);
            Cz_Librations(10*m+c+1,1) = PCtemp(1,919+30*m+c);
        }
    }
    if (0<=dt && dt<=8) {
        j=0;
        Mjd0 = t1;
    }
    else if(8<dt && dt<=16) {
        j=1;
        Mjd0 = t1+8*j;
    }
    else if (16<dt && dt<=24) {
        j=2;
        Mjd0 = t1+8*j;
    }
    else if(24<dt && dt<=32) {
        j=3;
        Mjd0 = t1+8*j;
    }
    double Cx_Librations2[10],Cy_Librations2[10],Cz_Librations2[10];
    for (int c = 0; c < 10; ++c) {
        Cx_Librations2[c] = Cx_Librations(j*10+c+1,1);
        Cy_Librations2[c] = Cy_Librations(j*10+c+1,1);
        Cz_Librations2[c] = Cz_Librations(j*10+c+1,1);
    }

    Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8,Cx_Librations2, Cy_Librations2,Cz_Librations2, res);
    Matrix Librations(3,1);
    for (int c = 0; c < 3; c++) {
        Librations(c+1,1) = res[c];
    }


    Matrix r_MoonAux(3,1);
    for (int c = 0; c < 3; c++) {
        r_MoonAux(c+1,1) = r_Moon(c+1,1)*EMRAT1;
    }

    Matrix r_EarthAux (3,1);

    r_Earth = r_Earth-r_MoonAux;
    for (int c = 0; c < 3; c++) {
        r_EarthAux(c+1,1)=r_Earth(c+1,1)*-1;
    }

    r_Mercury = r_EarthAux+r_Mercury;
    r_Venus = r_EarthAux+r_Venus;
    r_Mars = r_EarthAux+r_Mars;
    r_Jupiter = r_EarthAux+r_Jupiter;
    r_Saturn = r_EarthAux+r_Saturn;
    r_Uranus = r_EarthAux+r_Uranus;
    r_Neptune = r_EarthAux+r_Neptune;
    r_Pluto = r_EarthAux+r_Pluto;
    r_Sun = r_EarthAux+r_Sun;
}
