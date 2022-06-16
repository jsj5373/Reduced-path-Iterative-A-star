# include <iostream>
# include <eigen3/Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

int main(){
MatrixXd A = MatrixXd(3,3);
for (int i = 1 ; i<3 ; i++){
    A.col(i) << i,2*i,3*i;
}

cout << A << endl;
MatrixXd B = MatrixXd(3,2);
B.col(0) << (A.col(1) - A.col(0))/5;
cout<<B<<endl;
}