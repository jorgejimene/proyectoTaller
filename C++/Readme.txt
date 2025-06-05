Instrucciones para compliación de tests unitarios:
g++ -I../INCLUDE ../EKF_Test.cpp DEInteg.cpp Accel.cpp AccelHarmonic.cpp AuxParam.cpp Cheb3D.cpp globalMatrix.cpp G_AccelHarmonic.cpp VarEqn.cpp Matrix.cpp Mjday.cpp R_x.cpp R_y.cpp R_z.cpp Legendre.cpp SAT_Const.cpp NutAngles.cpp Unit.cpp AzElPa.cpp Frac.cpp AccelPointMass.cpp IERS.cpp MeanObliquity.cpp TimeUpdate.cpp Mjday_TBD.cpp MeasUpdate.cpp Sign_.cpp NutMatrix.cpp PoleMatrix.cpp PrecMatrix.cpp timediff.cpp EccAnom.cpp elements.cpp EqnEquinox.cpp gmst.cpp gast.cpp GHAMatrix.cpp JPL_Eph_DE430.cpp LTC.cpp Vector.cpp -o ../bin/programa.exe
(Pueden sustituirse los nombres de los archivos por *.cpp)
Seguido de: 
cd ../bin
Y estando ya dentro de /bin:
programa.exe pause