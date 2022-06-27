#include <Integrator.hpp>

Eigen::MatrixXd f1(Eigen::MatrixXd x0) 

{

    Eigen::MatrixXd mat; 
    mat.resize(x0.rows(), x0.cols());
    mat(0, 0) = x0(0, 0)*2;
    mat(1, 0) = x0(1, 0)*2; 

    return mat; 

}

int main() 

{
    
    Eigen::MatrixXd b; 
    b.resize(2, 1); 
    b(0, 0) = 1.0;
    b(1, 0) = 2.0; 

    std::cout << b << "\n"; 

    std::cout << Integrator::integrate_step_rk4_autonomous<f1>(b, 0.01);

    return 0; 

}