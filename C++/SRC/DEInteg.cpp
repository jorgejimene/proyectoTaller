#include "../INCLUDE/DEInteg.h"

enum DE_STATE {
    DE_INIT = 1,
    DE_DONE = 2,
    DE_BADACC = 3,
    DE_NUMSTEPS = 4,
    DE_STIFF = 5,
    DE_INVPARAM = 6
};
using namespace std;
Matrix DEInteg(function<Matrix(double, Matrix)> f, double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y) {
    double eps = 2.22044604925031e-16;
    double twou = 2*eps;
    double fouru = 4*eps;

    DE_STATE State_ = DE_STATE::DE_INIT;
    bool PermitTOUT = true;
    double told = 0;

    array<double,14> two = {1.0,2.0,4.0,8.0,16.0,32.0,64.0,128.0,256.0,512.0,1024.0,2048.0,4096.0,8192.0};
    array<double, 14> gstr = {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188, 0.0143, 0.0114, 0.00936, 0.00789, 0.00679, 0.00592, 0.00524, 0.00468};

    Matrix yy(n_eqn, 1), wt(n_eqn, 1),p(n_eqn,1),yp(n_eqn,1),phi(n_eqn,17),g(14,1),sig(14,1),rho(14,1),w(13,1),alpha(13,1),beta(13,1),v(13,1),psi_(13,1);

    if (t==tout) {
        return y;
    }

    double epsilon = max(relerr, abserr);
    if (relerr < twou || abserr < twou ||epsilon <= 0.0 || static_cast<int>(State_) > static_cast<int>(DE_STATE::DE_INVPARAM) ||(State_ != DE_STATE::DE_INIT && t!= told)) {
        State_ = DE_STATE::DE_INVPARAM;
    }
    double del = tout-t, absdel = abs(del);
    double tend = t+100.0*del;
    if (!PermitTOUT) tend = tout;

}