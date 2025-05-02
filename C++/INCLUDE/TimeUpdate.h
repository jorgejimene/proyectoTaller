//
// Created by jojimene on 10/04/2025.
//

#ifndef C___TIMEUPDATE_H
#define C___TIMEUPDATE_H
#include <Matrix.h>

/**
 * 
 * @param P 
 * @param Phi 
 * @param Qdt 
 * @return 
 */
void TimeUpdate(Matrix *P, Matrix *Phi, Matrix* Qdt=nullptr);
#endif //C___TIMEUPDATE_H
