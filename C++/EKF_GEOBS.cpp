#include <AzElPa.h>
#include <iomanip>
#include <iostream>
#include <Matrix.h>
#include <Position.h>
#include <TimeUpdate.h>

#include "Accel.h"
#include "AuxParam.h"
#include "DEInteg.h"
#include "globalMatrix.h"
#include "LTC.h"
#include "MeasUpdate.h"
#include "VarEqn.h"
/*
 * Archivo main
 */
int main() {
	Matrix eop = Matrix::LoadFromFile("eop19620101.txt");
    Matrix eopdata = eop.transpose();
    GGM03S();
    GEOS3();

    double sigma_range = 92.5;
    double sigma_az = 0.0224 * consts.Rad;
    double sigma_el = 0.0139 * consts.Rad;

    double lat = consts.Rad * 21.5748;
    double lon = consts.Rad * (-158.2706);
    double alt = 300.20;
    double rVec[3];
    Position(lon, lat, alt, rVec);
    Matrix pos(3, 1);
    pos(1,1) = rVec[0];
    pos(2,1) = rVec[1];
    pos(3,1) = rVec[2];
    if (pos.getFil() != 3 || pos.getCol() != 1) {
        std::cerr << R"(Error: pos tiene dimensiones incorrectas ()" << pos.getFil() << "x" << pos.getCol() << ")\n";
        return 1;
    }
    Matrix Rs = pos;
    if (Rs.getFil() != 3 || Rs.getCol() != 1) {
        std::cerr << "Error: Rs tiene dimensiones incorrectas (" << Rs.getFil() << "x" << Rs.getCol() << ")\n";
        return 1;
    }

    double Mjd1 = obs(1, 1);
    double Mjd2 = obs(9, 1);
    double Mjd3 = obs(18, 1);

    Matrix r2(3, 1), v2(3, 1);
    r2(1, 1) = 6221397.62857869;
    r2(2, 1) = 2867713.77965738;
    r2(3, 1) = 3006155.98509949;
    v2(1, 1) = 4645.04725161806;
    v2(2, 1) = -2752.21591588204;
    v2(3, 1) = -7507.99940987031;

    Matrix Y0_apr(6, 1);
    for (int i = 1; i <= 3; ++i) {
        Y0_apr(i, 1) = r2(i, 1);
        Y0_apr(i + 3, 1) = v2(i, 1);
    }

    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);

    double Mjd_UTC = obs(9, 1);
    auxParam.Mjd_UTC = Mjd_UTC;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.sun = true;
    auxParam.moon = true;
    auxParam.planets = true;
    int n_eqn = 6;

    Matrix Y = DEInteg(Accel, 0, -(Mjd_UTC - Mjd0) * 86400.0, 1e-13, 1e-6, 6, Y0_apr);
    if (Y.getFil() != 6 || Y.getCol() != 1) {
        std::cerr << "Error: Y tiene dimensiones incorrectas después de la primera integración (" << Y.getFil() << "x" << Y.getCol() << ")\n";
        return 1;
    }

    Matrix P(6, 6);
    for (int i = 1; i <= 3; ++i) P(i, i) = 1e8;
    for (int i = 4; i <= 6; ++i) P(i, i) = 1e3;

    Matrix LT = LTC(lon, lat);
    if (LT.getFil() != 3 || LT.getCol()!= 3) {
        std::cerr << "Error: LT tiene dimensiones incorrectas (" << LT.getFil() << "x" << LT.getCol() << ")\n";
        return 1;
    }

    Matrix yPhi(42, 1);
    Matrix Phi(6, 6);

    double t = 0;

    for (int i = 1; i <= 46; ++i) {
        double t_old = t;
        Matrix Y_old = Y;

        Mjd_UTC = obs(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;

        double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
        IERS(eopdata, Mjd_UTC,x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC,'l');
        double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
        timediff(UT1_UTC, TAI_UTC,UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

        double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
        double Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        auxParam.Mjd_UTC = Mjd_UTC;
        auxParam.Mjd_TT = Mjd_TT;

        for (int ii = 1; ii <= 6; ++ii) {
            yPhi(ii, 1) = Y_old(ii, 1);
            for (int j = 1; j <= 6; ++j) {
                yPhi(6 * j + ii, 1) = (ii == j) ? 1.0 : 0.0;
            }
        }

        yPhi = DEInteg(VarEqn, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);
        if (yPhi.getFil() != 42 || yPhi.getCol() != 1) {
            std::cerr << "Error: yPhi tiene dimensiones incorrectas en la iteración " << i << " (" << yPhi.getFil() << "x" << yPhi.getCol() << ")\n";
            return 1;
        }

        for (int j = 1; j <= 6; ++j) {
            for (int k = 1; k <= 6; ++k) {
                double val = yPhi(6 * j + k, 1);
                Phi(k, j) = val;
            }
        }

        Y = DEInteg(Accel, 0, t - t_old, 1e-13, 1e-6, 6, Y_old);
        if (Y.getFil() != 6 || Y.getCol()!= 1) {
            std::cerr << "Error: Y tiene dimensiones incorrectas en la iteración " << i << " (" << Y.getFil() << "x" << Y.getCol() << ")\n";
            return 1;
        }

        double theta = gmst(Mjd_UT1);
        Matrix U = R_z(theta);
        if (U.getFil() != 3 || U.getCol()!= 3) {
            std::cerr << "Error: U tiene dimensiones incorrectas (" << U.getFil() << "x" << U.getCol() << ")\n";
            return 1;
        }
        Matrix r(3, 1);
        for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
        Matrix Ur = U * r;
        if (Ur.getFil() != 3 || Ur.getCol() != 1) {
            std::cerr << "Error: U*r tiene dimensiones incorrectas (" << Ur.getFil() << "x" << Ur.getCol() << ")\n";
            return 1;
        }
        if (Ur.getCol() != Rs.getFil() || Ur.getCol() != Rs.getFil()) {
            std::cerr << "Error: Dimensiones incompatibles para U*r - Rs en la iteración " << i
                      << " (U*r: " << Ur.getFil() << "x" << Ur.getCol() << ", Rs: " << Rs.getFil() << "x" << Rs.getCol() << ")\n";
            return 1;
        }
        Matrix s = LT * (Ur - Rs);
        if (s.getFil() != 3 || s.getCol() != 1) {
            std::cerr << "Error: s tiene dimensiones incorrectas (" << s.getFil() << "x" << s.getCol() << ")\n";
            return 1;
        }

        TimeUpdate(&P, &Phi);
        double sV[3];
        sV[0] =s(1,1);
        sV[1] =s(2,1);
        sV[2] =s(3,1);
        double Azim, Elev;
        double dAds[3], dEds[3];
        Matrix dAdY(1,6);
        AzElPa(sV,Azim, Elev, dAds, dEds);
        Matrix new_dAdY(dAdY.getFil(), dAdY.getCol() + 3);

        for (int i = 1; i <= dAdY.getFil(); i++) {
            for (int j = 1; j <= dAdY.getCol(); j++) {
                new_dAdY(i, j) = dAdY(i, j);
            }
        }
        dAdY = new_dAdY;
        Matrix K(6,1), Y_new(6,1), P_new(6,6);
        Matrix z(1, 1), g(1, 1), s_sigma(1, 1);
        z(1, 1) = obs(i, 2);
        g(1, 1) = Azim;
        s_sigma(1, 1) = sigma_az;
        MeasUpdate(Y, z, g, s_sigma, dAdY, P, 6,K);
        Y = Y_new;
        P = P_new;

        for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
        Ur = U * r;
        s = LT * (Ur - Rs);

        AzElPa(sV, Azim, Elev,dAds, dEds);
        Matrix result(LT.getFil(), LT.getCol());
        for (int i=0;i<3;i++) {
            LT(i,1) = LT(i,1)*dEds[i];
        }
        Matrix dEdY = result * U;

        z(1, 1) = obs(i, 3);
        g(1, 1) = Elev;
        s_sigma(1, 1) = sigma_el;
        MeasUpdate(Y, z, g, s_sigma, dEdY, P, 6,K);
        Y = Y_new;
        P = P_new;

        for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
        Ur = U * r;
        s = LT * (Ur - Rs);
        double Dist = s.norm();
        for (int i=1;i<s.getCol();i++) {
            for (int j=1;j<s.getFil();i++) {
                s(i,j) = s(i,j) / Dist;
            }
        }
        Matrix dDds = s.transpose();
        Matrix dDdY = dDds * LT * U;
        z(1, 1) = obs(i, 4);
        g(1, 1) = Dist;
        s_sigma(1, 1) = sigma_range;
        MeasUpdate(Y, z, g, s_sigma, dDdY, P, 6,K);
        Y = Y_new;
        P = P_new;
    }

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    IERS(eopdata, obs(46, 1), x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC,'l');
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    auxParam.Mjd_UTC = Mjd_UTC;
    auxParam.Mjd_TT = Mjd_TT;

    Matrix Y0 = DEInteg(Accel, 0, -(obs(46, 1) - obs(1, 1)) * 86400.0, 1e-13, 1e-6, 6, Y);

    double aux[] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};
    Matrix Y_true(6, 1);
    for (int i = 0; i < 6; ++i) {
        Y_true(i + 1, 1) = aux[i];
    }

    std::cout << "\nError de Estimación de Posición\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "dX " << std::setw(10) << Y0(1, 1) - Y_true(1, 1) << " [m]\n";
    std::cout << "dY " << std::setw(10) << Y0(2, 1) - Y_true(2, 1) << " [m]\n";
    std::cout << "dZ " << std::setw(10) << Y0(3, 1) - Y_true(3, 1) << " [m]\n";
    std::cout << "\nError de Estimación de Velocidad\n";
    std::cout << "dVx " << std::setw(8) << Y0(4, 1) - Y_true(4, 1) << " [m/s]\n";
    std::cout << "dVy " << std::setw(8) << Y0(5, 1) - Y_true(5, 1) << " [m/s]\n";
    std::cout << "dVz " << std::setw(8) << Y0(6, 1) - Y_true(6, 1) << " [m/s]\n";

    return 0;
}