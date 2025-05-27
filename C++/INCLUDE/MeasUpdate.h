//
// Created by Jorge on 02/05/2025.
//

#ifndef MEASUPDATE_H
#define MEASUPDATE_H
#include <Matrix.h>

void MeasUpdate(Matrix &x, Matrix &z, Matrix &g, Matrix &s,
                Matrix &G, Matrix &P, int n, Matrix &K);
#endif //MEASUPDATE_H
