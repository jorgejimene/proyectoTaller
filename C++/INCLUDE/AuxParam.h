//
// Created by Jorge on 17/05/2025.
//

#ifndef AUXPARAM_H
#define AUXPARAM_H
struct AuxParam {
   mutable double Mjd_UTC;
   mutable int n;
   mutable int m;
   mutable int sun;
   mutable int moon;
   mutable int planets;
   mutable double Mjd_TT;
};
extern const AuxParam auxParam;
#endif //AUXPARAM_H
