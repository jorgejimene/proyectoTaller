//
// Created by Jorge on 17/05/2025.
//

#ifndef AUXPARAM_H
#define AUXPARAM_H
struct AuxParam {
   double Mjd_UTC;
   int n;
   int m;
   int sun;
   int moon;
   int planets;
   double Mjd_TT;
};
extern const AuxParam auxParam;
#endif //AUXPARAM_H
