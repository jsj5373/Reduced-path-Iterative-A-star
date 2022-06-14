#ifndef TrackCtrl_H
#define TrackCtrl_H
#include <iostream>
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
class TrackCtrl{
private:

    VectorXd x = VectorXd(3);
    VectorXd dx = VectorXd(3);
    VectorXd ie = VectorXd(2);
    VectorXd h_d = VectorXd(2);
    VectorXd dh_d = VectorXd(2);
    VectorXd ddh_d = VectorXd(2);
    double dt;
    MatrixXd M = MatrixXd(3,3);
    const char *mode;
    VectorXd constant = VectorXd(2);



public:
    VectorXd u = VectorXd(2);
    VectorXd e = VectorXd(2);
    VectorXd dh = VectorXd(2);
    VectorXd de = VectorXd(2);

    VectorXd Delta_par = VectorXd(3);
    MatrixXd Omega_per = MatrixXd(2,3);
    MatrixXd Delta_per;
    MatrixXd Omega_par;
    MatrixXd Delta = MatrixXd(3,3);
    MatrixXd dDelta = MatrixXd(3,3);
    MatrixXd R = MatrixXd(2,2);
    VectorXd h = VectorXd(2); //fixed
    MatrixXd D_per = MatrixXd(3,2);
    VectorXd uh = VectorXd(2);

    double l=0.5;
    double den1;
    double den2;
    double den3;
    double den4;
    double den5;
    double delt;
    double ddelt;
    MatrixXd trackctrl(VectorXd x, VectorXd dx, VectorXd ie, VectorXd h_d, VectorXd dh_d, VectorXd ddh_d,double dt,MatrixXd M, const char *mode, VectorXd constant){
            MatrixXd result(2,4);
    double th=x[2]; // theta
    double M1=M(0,0);
    double M2=M(1,1);
    double M3=M(2,2);

    MatrixXd Mdecomp(3,3);
    MatrixXd Cdecomp(3,3);
    
    Delta_par << l*sin(th), -l*cos(th), 1;
    Delta_par = Delta_par/(pow(l,2)+1);

    Omega_per << 1, 0, -l*sin(th), 0, 1, l*cos(th);
    
    Delta_per = (M.inverse()*Omega_per.transpose())*((Omega_per*M.inverse()*Omega_per.transpose()).inverse());

    Omega_par = (Delta_par.transpose()*M*Delta_par).inverse()*Delta_par.transpose()*M;
    // den1 = pow(M2+4*M3+(M1-M2)*pow(sin(th),2),2);
    // den2 = pow(4*M3+M2*pow(cos(th),2)+M1*pow(sin(th),2),2);
    den1 = 2*M2*pow(l,2)*cos(th)*sin(th);
    den2 = 2*M1*pow(l,2)*cos(th)*sin(th);//den4, den5 order fixed : should be declared before den3
    den4 = M1*pow(l,2)*pow(sin(th),2);
    den5 = M2*pow(l,2)*pow(cos(th),2);
    den3 = pow(den5+den4+M3,2);

    Delta << Delta_par, Delta_per;
    dDelta << l*cos(th)/(pow(l,2)+1), -(den5 + M3)*(den2 - den1)/den3 - den1/(den5 + den4 + M3), (den5 - M2*pow(l,2)*pow(sin(th),2))/(den5 + den4 + M3) - M2*pow(l,2)*cos(th)*sin(th)*(den2-den1)/den3,
                l*sin(th)/(pow(l,2)+1), (M1*pow(l,2)*pow(cos(th),2)-den4)/(den5+den4+M3) - M1*pow(l,2)*cos(th)*sin(th)*(den2-den1)/den3, den2/(den5+den4+M3) - (den4+M3)*(den2-den1)/den3,
                0,                      M1*l*sin(th)*(den2-den1)/den3 - M1*l*cos(th)/(den5+den4+M3), -M2*l*sin(th)/(den5+den4+M3) - M2*l*cos(th)*(den2-den1)/den3;
    
    dDelta = dDelta*dx(2);
    //cout<<dDelta<<endl;
    D_per = Delta.block(0,1,3,2);
    R << cos(th), -sin(th), sin(th), cos(th);
    h << x[0]+l*cos(th), x[1]+l*sin(th); 
    dh = Omega_per*dx;
    e = h - h_d;
    de = dh-dh_d;
    ie = ie + dt*e;

    if ((strcmp(mode, "pidconst"))==0){
        double K_p = constant(0);
        double K_i = constant(1);
        double K_v = constant(2);
        uh = -K_v*de-K_p*e-K_i*ie;
    }
    else{
        Mdecomp = Delta.transpose()*M*Delta;
        Cdecomp = Delta.transpose()*(M*dDelta);
        // cout<<"mdcmp"<<Mdecomp<<endl;
        // cout<<"cdcmp"<<Cdecomp<<endl;
        // cout<<"e"<<e<<endl;
        // cout<<"de"<<de<<endl;
        // cout<<"ie"<<ie<<endl;
        // cout<<"dh"<<dh<<endl;
        
        if (strcmp(mode,"pid")==0){
            double lambda = constant(0);
            MatrixXd K_v = 3*lambda*Mdecomp.block(1,1,2,2)-Cdecomp.block(1,1,2,2);
            MatrixXd K_p = 3*pow(lambda,2)*Mdecomp.block(1,1,2,2);
            MatrixXd K_i = pow(lambda,3)*Mdecomp.block(1,1,2,2);
            uh = -K_v*de-K_p*e-K_i*ie + Mdecomp.block(1,1,2,2)*ddh_d+Cdecomp.block(1,1,2,2)*dh_d;
        }
        else if (strcmp(mode,"backstepping")==0){
            double k1 = constant(0);
            double k2 = constant(1);
            uh = Cdecomp.block(1,1,2,2)*dh+Mdecomp.block(1,1,2,2)*(ddh_d-(1+k1*k2)*e-(k1+k2)*de);
            //cout<<"uh"<<uh<<endl;
        }
    }
    Eigen::MatrixXd tmp(3,2);
    tmp << cos(x(2)), 0, sin(x(2)), 0, 0, 1;
    u = (D_per.transpose()*tmp).inverse()*uh;
    result << u, ie, e, dh;
    return result;
    }

};

#endif