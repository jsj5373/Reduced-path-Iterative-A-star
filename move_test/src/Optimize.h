# ifndef Optimize_H
# define Optimize_H
# include <iostream>
# include <eigen3/Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
class Optimize {
private:

    
public: 
    double dt;
    MatrixXd points; 
    VectorXd V_init = VectorXd(3); 
    VectorXd A_init = VectorXd(3); 

    int q_size = points.cols()-1; // 총 분할 구간 수
    MatrixXd q = points;
    

    int m = 10; //구간당 분할 수
    double dT = m*dt; //구간당 소요시간
    int N =  q_size * m;
    
    MatrixXd q_V = MatrixXd(3,q_size+1);
    MatrixXd q_A = MatrixXd(3,q_size+1);

    

    MatrixXd delta_x = MatrixXd(3,q_size);
    MatrixXd delta_y = MatrixXd(3,q_size);
    MatrixXd delta_z = MatrixXd(3,q_size);
    MatrixXd coeff_x = MatrixXd(3,q_size+1);
    MatrixXd coeff_y = MatrixXd(3,q_size+1);
    MatrixXd coeff_z = MatrixXd(3,q_size+1);

    MatrixXd s_t = MatrixXd(3,N+1);
    MatrixXd v_t = MatrixXd(3,N+1);
    MatrixXd a_t = MatrixXd(3,N+1);
    
    s_t.block(1,1,3,1) = q.block(1,1,3,1); //set start point 
    
 



};


int main(){


    for (int i = 1; i < q_size+2 ; i++)
        
        if (i == 1)
            q_V.block(1,i,3,i) = V_init; 
            q_A.block(1,i,3,i) = A_init; //이 eigen 구문 대입 되는지 check
        else
            q_V.block(1,i,3,i) = 2*(q.block(1,i,3,i)-q.block(1,i-1,3,i-1))/dT-q_V.block(1,i-1,3,i-1); 
            q_A.block(1,i,3,i) = 2*(q_V.block(1,i,3,i)-q_V.block(1,i-1,3,i-1)/dT)-q_A.block(1,i-1,3,i-1); // 이것도 check
        
    for (i = 1; i < q_size+1 ;i++)
        delta_x.block(1,i,3,i) << q(1,i+1)-q(1,i)-q_V(1,i)*dT-(1/2)*q_A(1,i)*(dT^2), q_V(1,i+1)-q_V(1,i)-q_A(1,i)*dT, q_A(1,i+1)-q_A(1,i);
        delta_y.block(1,i,3,i) << q(2,i+1)-q(2,i)-q_V(2,i)*dT-(1/2)*q_A(2,i)*(dT^2), q_V(2,i+1)-q_V(2,i)-q_A(2,i)*dT, q_A(2,i+1)-q_A(2,i);
        delta_z.block(1,i,3,i) << q(3,i+1)-q(3,i)-q_V(3,i)*dT-(1/2)*q_A(3,i)*(dT^2), q_V(3,i+1)-q_V(3,i)-q_A(3,i)*dT, q_A(3,i+1)-q_A(3,i);
        MatrixXd Temp = MatrixXd(3,3);
        Temp << 720, -360*dT, 60*pow(dT,2), -360*dT, 168*pow(dT,2), -24*pow(dT,3), 60*pow(dT,2), -24*pow(dT,3), 3*pow(dT,4); 
        coeff_x.block(1,i,3,i) = (1/pow(dT,5))*Temp*delta_x.block(1,i,3,i);
        coeff_y.block(1,i,3,i) = (1/pow(dT,5))*Temp*delta_y.block(1,i,3,i);
        coeff_z.block(1,i,3,i) = (1/pow(dT,5))*Temp*delta_z.block(1,i,3,i);





    for (int i = 1 ; i <= N+1 ; ++i )
        t = ((i-1) % m) * dt;
        int j;
        j = 1 + (i-1)/m ;

        s_t(1,i) = pow(t,5)*coeff_x(1,j)/120 + pow(t,4)*coeff_x(2,j)/24 + pow(t,3)*coeff_x(3,j)/6 + (1/2)*q_A(1,j)*pow(t,2) + q_V(1,j)*t + q(1,j);
        s_t(2,i) = pow(t,5)*coeff_y(1,j)/120 + pow(t,4)*coeff_y(2,j)/24 + pow(t,3)*coeff_y(3,j)/6 + (1/2)*q_A(2,j)*pow(t,2) + q_V(2,j)*t + q(2,j);
        s_t(3,i) = pow(t,5)*coeff_z(1,j)/120 + pow(t,4)*coeff_z(2,j)/24 + pow(t,3)*coeff_z(3,j)/6 + (1/2)*q_A(3,j)*pow(t,2) + q_V(3,j)*t + q(3,j); 
        v_t(1,i) = pow(t,4)*coeff_x(1,j)/24 + pow(t,3)*coeff_x(2,j)/6 + pow(t,2)*coeff_x(3,j)/2 + q_A(1,j)*t + q_V(1,j);  
        v_t(2,i) = pow(t,4)*coeff_y(1,j)/24 + pow(t,3)*coeff_y(2,j)/6 + pow(t,2)*coeff_y(3,j)/2 + q_A(2,j)*t + q_V(2,j);
        v_t(3,i) = pow(t,4)*coeff_z(1,j)/24 + pow(t,3)*coeff_z(2,j)/6 + pow(t,2)*coeff_z(3,j)/2 + q_A(3,j)*t + q_V(3,j);
        a_t(1,i) = pow(t,3)*coeff_x(1,j)/6 + pow(t,2)*coeff_x(2,j)/2 + t*coeff_x(3,j) + q_A(1,j);
        a_t(2,i) = pow(t,3)*coeff_y(1,j)/6 + pow(t,2)*coeff_y(2,j)/2 + t*coeff_y(3,j) + q_A(2,j);
        a_t(3,i) = pow(t,3)*coeff_z(1,j)/6 + pow(t,2)*coeff_z(2,j)/2 + t*coeff_z(3,j) + q_A(3,j);
     
}



#endif
