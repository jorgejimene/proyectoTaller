#include "../INCLUDE/DEInteg.h"


/*
%----------------------------------------------------------------------------
%
% Purpose:
%   Numerical integration methods for ordinaray differential equations
%
%   This module provides implemenation of the variable order variable
%   stepsize multistep method of Shampine & Gordon.
%
% Last modified:   2015/08/25   M. Mahooti
%
% Reference:
%
%   Shampine, Gordon: "Computer solution of Ordinary Differential Equations",
%   Freeman and Comp., San Francisco (1975).
%
%----------------------------------------------------------------------------
*/

using namespace std;


Matrix DEInteg(Matrix (*func)(double x, const Matrix &Y), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y) {
    if (y.getCol() >1) {
        y = y.transpose();
    }

    int contador=0;
    double twou= 0.,fouru= 0.,epsilon= 0.,del= 0.,absdel= 0.,tend= 0.,releps= 0.,abseps= 0.,x= 0.,delsgn= 0.,h= 0.,hi = 0.,
        temp1= 0.,term= 0.,psijm1= 0.,gamma= 0.,eta= 0.,p5eps= 0.,round= 0.,sum= 0.,absh= 0.,hold= 0.,hnew= 0.,temp2= 0.,
        temp4= 0.,temp5= 0.,temp3= 0.,temp6= 0.,tau= 0.,xold= 0.,erkm2= 0.,erkm1= 0.,erk= 0.,err= 0.,auxrho= 0.,erkp1= 0.,
        r= 0.,told= 0.;
    int State_= 0,nostep= 0,kle4= 0,kold = 0,ki = 0,ifail = 0,k = 0,kp1 = 0,kp2 = 0,km1 = 0,km2 = 0,ns = 0,nsp1 = 0,
        realns = 0,im1 = 0,reali = 0,nsm2 = 0,limit1 = 0,nsp2 = 0,limit2 = 0,ip1 = 0,knew = 0;
    bool PermitTOUT = false,stiff = false,start = false,OldPermit = false,crash = false, phase1 = false,nornd = false,success = false;

    struct{
        int DE_INIT = 1;
        int DE_DONE = 2;
        int DE_BADACC = 3;
        int DE_NUMSTEPS = 4;
        int DE_STIFF = 5;
        int DE_INVPARAM = 6;
    }DE_STATE;

    twou = 2.0 * 2.2204e-16;
    fouru = 4.0 * 2.2204e-16;

    State_ = DE_STATE.DE_INIT;
    PermitTOUT = true;
    told = 0.0;

    double twoVec[] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0};
    double gstrVec[] = {1.0, 0.5, 0.0833, 0.0417, 0.0264, 0.0188, 0.0143, 0.0114, 0.00936, 0.00789, 0.00679, 0.00592, 0.00524, 0.00468};
    Matrix two(1,14,twoVec,14);
    Matrix gstr(1,14,gstrVec,14);
    Matrix yy(n_eqn, 1), wt(n_eqn, 1), p(n_eqn, 1), yp(n_eqn, 1), phi(n_eqn, 17);
    Matrix g(14, 1), sig(14, 1), rho(14, 1), w(13, 1), alpha(13, 1);
    Matrix beta(13, 1), v(13, 1), psi_(13, 1);

    if (t == tout) {
        return y;
    }

    epsilon = max(relerr, abserr);

    if (relerr < twou || abserr < twou || epsilon <= 0.0 ||
        State_ > DE_STATE.DE_INVPARAM ||
        (State_ != DE_STATE.DE_INIT && t != told)) {
        State_ = DE_STATE.DE_INVPARAM;
        return y;
        }
    del = tout - t, absdel = fabs(del);
    tend = t + 100.0 * del;
    if (!PermitTOUT) {
        tend = tout;
    }

    nostep = 0;
    kle4 = 0;
    stiff = false;
    releps = relerr/epsilon;
    abseps = abserr/epsilon;

    if (State_ == DE_STATE.DE_INIT || (!OldPermit) || (delsgn*del <=0.0)){
        start = true;
        x = t;
        yy = y;
        delsgn = sign_(1.0,del);
        h = sign_(max(fouru*abs(x),abs(tout - x)), tout - x);
    }

    while (true) {
        cout << "iteracion" << contador << endl;
        contador++;
        /*cout << "fabs(x-t)" << fabs(x-t) << endl;
        cout << "absdel" << absdel << endl;
        */
        if (fabs(x - t) >= absdel) {
            Matrix yout(n_eqn, 1), ypout(n_eqn, 1);
            g(2, 1) = 1.0;
            rho(2, 1) = 1.0;
            hi = tout - x;
            ki = kold + 1;

            for (int i = 0; i <= ki; ++i) {
                temp1 = i;
                w(i + 1, 1) = 1.0 / i;
            }
            term = 0.0;
            for (int j = 2; j <= ki; ++j) {
                psijm1 = psi_(j, 1);
                gamma = (hi + term) / psijm1, eta = hi / psijm1;
                for (int i = 1; i <= ki+1-j; ++i)
                    w(i + 1, 1) = gamma * w(i + 1, 1) - eta * w(i + 2, 1);
                g(j + 1, 1) = w(2, 1);
                rho(j + 1, 1) = gamma * rho(j, 1);
                term = psijm1;
            }

            for (int j = 1; j <= ki; ++j) {
                int i = ki - j;
                for (int l = 1; l <= n_eqn; ++l) {
                    yout(l, 1) = yout(l,1) + g(i + 1, 1) * phi(l, i + 1);
                    ypout(l, 1) = ypout(l,1) + rho(i + 1, 1) * phi(l, i + 1);
                }
            }
            yout = y + yout * hi;
            y = yout;
            State_ = DE_STATE.DE_DONE;
            t = tout;
            told = t;
            OldPermit = PermitTOUT;
            return y;
        }

        if (!PermitTOUT && (fabs(tout - x) < fouru * fabs(x))) {
            h = tout - x;
            yp = (*func)(x, yy);
            y = yy + yp * h;
            State_ = DE_STATE.DE_DONE;
            t = tout;
            told = t;
            OldPermit = PermitTOUT;
            return y;
        }

        h = sign_(min(fabs(h), fabs(tend - x)), h);
        for (int l = 1; l <= n_eqn; ++l) {
            wt(l, 1) = releps * fabs(yy(l, 1)) + abseps;
        }
        if (fabs(h) < fouru * fabs(x)) {
            h = sign_(fouru*fabs(x),h);
            crash = true;
            return y;
        }

        p5eps = 0.5 * epsilon;
        crash = false;
        g(2, 1) = 1.0;
        g(3, 1) = 0.5;
        sig(2, 1) = 1.0;

        ifail = 0;
        round = 0.0;
        for (int l = 1; l <= n_eqn; ++l){
            round = round + (y(l, 1) * y(l, 1)) / (wt(l, 1) * wt(l, 1));
        }
        round = twou * sqrt(round);
        if (p5eps < round) {
            epsilon = 2.0*round*(1.0+fouru);
            crash = true;
            return y;
        }

        if (start) {
            yp = (*func)(x, y);
            //yp.print();
            sum = 0.0;
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, 2) = yp(l, 1);
                phi(l, 3) = 0.0;
                sum += (yp(l, 1) * yp(l, 1)) / (wt(l, 1) * wt(l, 1));
            }
            sum = sqrt(sum);
            absh = abs(h);
            if (epsilon < 16.0 * sum * h * h){
                absh = 0.25 * sqrt(epsilon / sum);
            }
            h = sign_(max(absh, fouru * fabs(x)), h);
            hold = 0.0;
            hnew = 0.0;
            k = 1;
            kold = 0;
            start = false;
            phase1 = true;
            nornd = true;
            if (p5eps <= 100.0 * round) {
                nornd = false;
                for (int l = 1; l <= n_eqn; ++l) {
                    phi(l, 16) = 0.0;
                }
            }
        }
        while (true) {
            kp1 = k + 1;
            kp2 = k + 2;
            km1 = k - 1;
            km2 = k - 2;

            if (h != hold) {
                ns = 0;
            }
            if (ns <= kold) {
                ns++;
            }
            nsp1 = ns + 1;

            if (k >= ns) {
                beta(ns+1, 1) = 1.0;
                realns = ns;
                alpha(ns+1, 1) = 1.0 / realns;
                temp1 = h * realns;
                sig(nsp1+1, 1) = 1.0;
                if (k >= nsp1) {
                    for (int i = nsp1; i <= k; ++i) {
                        im1 = i - 1;
                        temp2 = psi_(im1+1, 1);
                        psi_(im1+1, 1) = temp1;
                        beta(i+1, 1) = beta(im1+1, 1) * psi_(im1+1, 1) / temp2;
                        temp1 = temp2 + h;
                        alpha(i+1, 1) = h / temp1;
                        reali = i;
                        sig(i + 2, 1) = reali * alpha(i+1, 1) * sig(i+1, 1);
                    }
                }
                psi_(k+1, 1) = temp1;

                if (ns > 1) {
                    if (k > kold) {
                        temp4 = k * kp1;
                        v(k+1, 1) = 1.0 / temp4;
                        nsm2 = ns - 2;
                        for (int j = 1; j <= nsm2; ++j) {
                            int i = k - j;
                            v(i+1, 1) = v(i+1, 1) - (alpha(j + 2, 1) * v(i + 2, 1));
                        }
                    }
                    limit1 = kp1 - ns;
                    temp5 = alpha(ns+1, 1);
                    for (int iq = 1; iq <= limit1; ++iq) {
                        v(iq+1, 1) = v(iq+1, 1) - (temp5 * v(iq + 2, 1));
                        w(iq+1, 1) = v(iq+1, 1);
                    }
                    g(nsp1+1, 1) = w(2, 1);
                } else {
                    for (int iq = 1; iq <= k; ++iq) {
                        temp3 = iq * (iq + 1);
                        v(iq+1, 1) = 1.0 / temp3;
                        w(iq+1, 1) = v(iq+1, 1);
                    }
                }

                nsp2 = ns + 2;
                if (kp1 >= nsp2) {
                    for (int i = nsp2; i <= kp1; ++i) {
                        limit2 = kp2 - i;
                        temp6 = alpha(i, 1);
                        for (int iq = 1; iq <= limit2; ++iq) {
                            w(iq+1, 1) =w(iq+1,1) - (temp6 * w(iq + 2, 1));
                        }
                        g(i+1, 1) = w(2, 1);
                    }
                }
            }

            if (k>=nsp1) {
                for (int i = nsp1; i <= k; ++i) {
                    temp1 = beta(i+1, 1);
                    for (int l = 1; l <= n_eqn; ++l) {
                        phi(l, i+1) = temp1 * phi(l, l+1);
                    }
                }
            }
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, kp2+1) = phi(l, kp1+1);
                phi(l, kp1+1) = 0.0;
                p(l, 1) = 0.0;
            }
            for (int j = 1; j <= k; ++j) {
                int i = kp1 - j;
                ip1 = i + 1;
                temp2 = g(i+1, 1);
                for (int l = 1; l <= n_eqn; ++l) {
                    p(l, 1) = p(l,1) + (temp2 * phi(l, i+1));
                    phi(l, i+1) = phi(l,i+1) + phi(l, ip1+1);
                }
            }
            if (nornd) {
                /*cout << "y" << endl;
                y.print();
                cout << "p" << endl;
                p.print();
                cout << "h" << h << endl;
                */
                p = y + p * h;
                //p.print();
            } else {
                for (int l = 1; l <= n_eqn; ++l) {
                    tau = h * p(l, 1) - phi(l, 16);
                    p(l, 1) = y(l, 1) + tau;
                    phi(l, 17) = (p(l, 1) - y(l, 1)) - tau;
                }
            }
            xold = x;
            x = x+ h;
            absh = abs(h);
            yp = (*func)(x, p);

            erkm2 = 0.0;
            erkm1 = 0.0;
            erk = 0.0;
            for (int l = 1; l <= n_eqn; ++l) {
                temp3 = 1.0 / wt(l, 1);
                temp4 = yp(l, 1) - phi(l, 1+1);
                wt.print();
                yp.print();
                phi.print();
                cout << "temp3" << temp3 << endl;
                cout << "temp4" << temp4 << endl;
                if (km2 > 0) {
                    erkm2 = erkm2 + ((phi(l,km1+1)+temp4)*temp3)*((phi(l,km1+1)+temp4)*temp3);
                }
                if (km2 >= 0) {
                    erkm1 = erkm1 + ((phi(l,k+1)+temp4)*temp3)*((phi(l,k+1)+temp4)*temp3);
                }
                erk = erk + (temp4*temp3)*(temp4*temp3);
            }

            if (km2 > 0) {
                erkm2 = absh * sig(km1 + 1, 1) * gstr(1,km2+1) * sqrt(erkm2);
            }
            if (km2 >= 0) {
                erkm1 = absh * sig(k + 1, 1) * gstr(1,km1+1) * sqrt(erkm1);
            }

            temp5 = absh * sqrt(erk);
            err = temp5 * (g(k+1, 1) - g(kp1+1, 1));
            erk = temp5 * sig(kp1, 1) * gstr(1,k + 1);
            knew = k;

            if (km2 > 0){
                if (max(erkm1,erkm2) <= erk) {
                    knew = km1;
                }
            }
            if (km2 == 0){
                if (erkm1 <=0.5*erk) {
                    knew = km1;
                }
            }


            success = (err <= epsilon);

            if (!success) {
                phase1 = false;
                x = xold;
                for (int i = 1; i <= k; ++i) {
                    temp1 = 1.0 / beta(i+1, 1);
                    ip1 = i + 1;
                    for (int l = 1; l <= n_eqn; ++l)
                        phi(l, i+1) = temp1 * (phi(l, i+1) - phi(l, ip1+1));
                }
                if (k >= 2) {
                    for (int i = 2; i <= k; ++i) {
                        psi_(i, 1) = psi_(i + 1, 1) - h;
                    }
                }
                ifail= ifail+1;
                temp2 = 0.5;
                if (ifail > 3){
                    if (p5eps < (0.25 * erk)) {
                        temp2 = sqrt(p5eps / erk);
                    }
                }
                if (ifail >= 3) {
                    knew = 1;
                }
                h = temp2*h;
                k = knew;
                if (fabs(h) < fouru * fabs(x)) {
                    crash = true;
                    h = sign_(fouru * fabs(x),h);
                    epsilon = epsilon * 2.0;
                    return y;
                }
            }
            //p.print();
            if (success) {
                break;
            }

        }


        kold = k;
        hold = h;

        temp1 = h * g(kp1, 1);
        if (nornd) {
            for (int i=1; i<=n_eqn; ++i) {
                y(i,1) = p(i,1) + temp1 * (yp(i,1) - phi(i,2));
            }
        }
        else {
            for (int i=1; i<=n_eqn; ++i) {
                auxrho = (temp1 * (yp(i,1)-phi(i,2)) - phi(i,17));
                y(i,1) = p(i,1) + auxrho;
                phi(i,16) = (y(i,1) - p(i,1)) - auxrho;
            }
        }
        yp = (*func)(x, y);
        //yp.print();

        for (int l = 1; l <= n_eqn; ++l) {
            phi(l, kp1+1) = yp(l, 1) - phi(l, 2);
            phi(l, kp2+1) = phi(l, kp1+1) - phi(l, kp2+1);
        }
        for (int i = 1; i <= k; ++i) {
            for (int l = 1; l <= n_eqn; ++l) {
                phi(l, i+1) = phi(l,i+1) + phi(l, kp1+1);
            }
        }

        erkp1 = 0.0;
        if (knew == km1 || k == 12) phase1 = false;

        if (phase1) {
            k = kp1;
            erk = erkp1;
        } else {
            if (knew == km1) {
                k = km1;
                erk = erkm1;
            } else {
                if (kp1 <= ns) {
                    for (int l = 1; l <= n_eqn; ++l) {
                        erkp1 =  erkp1 + (phi(l,kp2+1)/wt(l,1))*(phi(l,kp2+1)/wt(l,1));
                    }
                    erkp1 = absh * gstr(1,kp1+1) * sqrt(erkp1);
                    if (k > 1) {
                        if (erkm1 <= min(erk, erkp1)) {
                            k = km1;
                            erk = erkm1;
                        } else {
                            if (erkp1 < erk && k != 12) {
                                k = kp1;
                                erk = erkp1;
                            }
                        }
                    } else if (erkp1 < 0.5 * erk) {
                        k = kp1;
                        erk = erkp1;
                    }
                }
            }
        }

        if (phase1 || p5eps >= erk * two(1,k+2)) {
            hnew = 2.0 * h;
        } else {
            if (p5eps < erk) {
                temp2 = k + 1;
                r = pow(p5eps / erk, 1.0 / temp2);
                hnew = absh * max(0.5, min(0.9, r));
                hnew = sign_(max(hnew, fouru * fabs(x)), h);
            } else {
                hnew = h;
            }
        }
        h = hnew;

        if (crash) {
            State_ = DE_STATE.DE_BADACC;
            relerr = epsilon*releps;
            abserr = epsilon*abseps;
            y = yy;
            t = x;
            told = t;
            OldPermit = true;
            return y;
        }

        nostep++;
        kle4++;
        if (kold > 4) {
            kle4 = 0;
        }
        if (kle4 >= 500) {
            stiff = true;
        }
    }
}
