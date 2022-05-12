#include "MathUtils.hpp"

using namespace MathUtils;
#define assertm(exp, msg) assert(((void)msg, exp))

void MathUtils::Timer::tic() 

{

    t_start = std::chrono::high_resolution_clock::now();

}

void MathUtils::Timer::toc() {

    t_end = std::chrono::high_resolution_clock::now();
    elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
    std::cout << "Time elapsed: "
              << elapsed_time_ms << " milliseconds" << std::endl;

}

Eigen::Matrix<double, 3, 1> MathUtils::so3toVec(Eigen::Matrix3d so3mat) {

    Eigen::Matrix<double, 3, 1> Vec; 
    Vec << so3mat(2,1), so3mat(0, 2), so3mat(1, 0);
    return Vec;

}

Eigen::Matrix<double, 6, 1> MathUtils::se3toVec(Eigen::Matrix4d se3mat) {

    Eigen::Matrix<double, 6, 1> Vec; 
    Vec <<  se3mat(0, 3), se3mat(1, 3), se3mat(2, 3), se3mat(2,1), se3mat(0,2), se3mat(1, 0);

    return Vec;

}

Eigen::MatrixXd MathUtils::forwardFiniteDifferences(Eigen::MatrixXd mat, Eigen::MatrixXd mat_plus, double eps) {

    Eigen::MatrixXd derivative;
    assertm(mat.rows() == mat_plus.rows(), "mat and mat_plus must be the same size"); 
    assertm(mat.cols() == mat_plus.cols(), "mat and mat_plus must be the same size");
    derivative.resize(mat.rows(), mat.cols());

    derivative = (mat_plus - mat).array() / eps; 

    std::cout.precision(5);

    // std::cout << "mat_plus: " << derivative << std::endl;
    // std::cout << "mat: " << mat << std::endl;

    return derivative; 

}

Eigen::Matrix<double, 3, 3> MathUtils::skew_v(Eigen::Vector<double, 3> v) 

{

    Eigen::Matrix<double, 3, 3> skew_mat;

    skew_mat << 0, - v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;

    return skew_mat; 
    

}

Eigen::Matrix<double, 3, 3> MathUtils::skew_m(Eigen::Matrix<double, 3, 1> v) 

{

    Eigen::Matrix<double, 3, 3> skew_mat;

   skew_mat << 0, -v(2, 0), v(1, 0), v(2, 0), 0, -v(0, 0), -v(1, 0), v(0, 0), 0;

    return skew_mat; 
    

}


Eigen::Matrix<double, 6, 6> MathUtils::adjointTransformation(Eigen::Affine3d T) 

{

    Eigen::Matrix<double, 6, 6> adjT; 

    adjT.setZero();
    adjT.block<3, 3>(0, 0) = T.rotation();
    adjT.block<3, 3>(3, 3) = T.rotation();
    adjT.block<3, 3>(0, 3) = skew_v(T.translation())*T.rotation();

    return adjT;

}
