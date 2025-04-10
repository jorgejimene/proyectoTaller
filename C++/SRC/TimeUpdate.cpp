//
// Created by jojimene on 10/04/2025.
//
#include "../INCLUDE/TimeUpdate.h"

double TimeUpdate(double P, double Phi, double Qdt){
    //hay que cambiar el segundo Phi por Phi'
    return Phi*P*Phi + Qdt;
}