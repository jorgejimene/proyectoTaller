#include "../INCLUDE/Matrix.h"
#include "../INCLUDE/SAT_Const.h"

/*
%--------------------------------------------------------------------------
%
% IERS: Management of IERS time and polar motion data
%
% Last modified:   2018/02/01   M. Mahooti
%
%--------------------------------------------------------------------------
*/
void IERS(Matrix &eop, double Mjd_UTC,  double &x_pole, double &y_pole, double &UT1_UTC, double &LOD, double &dpsi, double &deps, double &dx_pole, double &dy_pole, double &TAI_UTC,char interp = 'o'){


    if(interp == 'o') {
        interp = 'n';
    }

    if (interp == 'l') {
        // linear interpolation
        int mjd = (int)(Mjd_UTC);
        int i = -1;
        int filas = eop.getFil();
        int columnas = eop.getCol();
        for (int j = 1; j< columnas; j++){
            if(mjd == (int) eop(4,j)){
                i = j;
                break;
            }
        }
        double mfme = 1440 * (Mjd_UTC - (int)(Mjd_UTC));
        double fixf = mfme / 1440;
        //Setting of IERS Earth rotation parameters
        // (UT1 - UTC[s], TAI - UTC[s], x["], y ["])
        x_pole = eop(5,i) + (eop(5,i+1) - eop(5,i)) * fixf;
        y_pole = eop(6,i) + (eop(6,i+1) - eop(6,i)) * fixf;
        UT1_UTC = eop(7,i) + (eop(7,i+1) - eop(7,i)) * fixf;
        LOD = eop(8,i) + (eop(8,i+1) - eop(8,i)) * fixf;
        dpsi = eop(9,i) + (eop(9,i+1) - eop(9,i)) * fixf;
        deps = eop(10,i) + (eop(10,i+1) - eop(10,i)) * fixf;
        dx_pole = eop(11,i) + (eop(11,i+1) - eop(11,i)) * fixf;
        dy_pole = eop(12,i) + (eop(12,i+1) - eop(12,i)) * fixf;
        TAI_UTC = eop(13,i);

        x_pole = x_pole /consts.Arcs;  // Pole coordinate[rad]
        y_pole = y_pole /consts.Arcs;  // Pole coordinate[rad]
        dpsi = dpsi /consts.Arcs;
        deps = deps /consts.Arcs;
        dx_pole = dx_pole /consts.Arcs; // Pole coordinate[rad]
        dy_pole = dy_pole /consts.Arcs; // Pole coordinate[rad]
    }
    else {
            int mjd = (int)(Mjd_UTC);
            int i = -1;
            int filas = eop.getFil();
            int columnas = eop.getCol();
            for (int j = 1; j< columnas; j++){
                if(mjd == (int) eop(4,j)){
                    i = j;
                    break;
                }
            // Setting of IERS Earth rotation parameters
            // (UT1 - UTC[s], TAI - UTC[s], x["], y ["])
            x_pole = eop(5,i) /consts.Arcs;  // Pole coordinate[rad]
            y_pole = eop(6,i) /consts.Arcs;  // Pole coordinate[rad]
            UT1_UTC = eop(7,i);             // UT1 - UTC time difference[s]
            LOD = eop(8,i);             // Length of day[s]
            dpsi = eop(9,i) /consts.Arcs;
            deps = eop(10,i) /consts.Arcs;
            dx_pole = eop(11,i) /consts.Arcs; // Pole coordinate[rad]
            dy_pole = eop(12,i) /consts.Arcs; // Pole coordinate[rad]
            TAI_UTC = eop(13,i);            // TAI - UTC time difference[s]
        }
    }
}
