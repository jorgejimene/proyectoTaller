//
// Created by jojimene on 10/04/2025.
//
#include "../INCLUDE/TimeUpdate.h"

void TimeUpdate(Matrix *P, Matrix *Phi,Matrix* Qdt){
    Matrix Phi_t = Phi->transpose();

    *P = (*Phi) * (*P) * Phi_t;

    if (Qdt!=nullptr) {
        *P = *P + (*Qdt);
    }
}