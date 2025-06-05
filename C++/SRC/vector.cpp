//
// Created by jojimene on 08/05/2025.
//

#include "../INCLUDE/vector.h"
double norm(const double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
double dot(const double r[], const double x[]){
    return(r[0]*x[0]+r[1]*x[1]+r[2]*x[2]);
}
void cross(const double v1[], const double v2[], double v[]){
    v[0] = v1[1]* v2[2] - v1[2] * v2[1];
    v[1] = v1[2] * v2[0] - v1[0] *v2[2];
    v[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
