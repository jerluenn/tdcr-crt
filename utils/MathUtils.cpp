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

Eigen::Quaterniond MathUtils::rot2quat(Eigen::Matrix3d rot) 

{

    Eigen::Quaterniond eta(rot);

    return eta; 

}

Eigen::Matrix3d MathUtils::quat2Rot(Eigen::Matrix<double, 4, 1> eta) 

{

    Eigen::Matrix3d R;

    R(0,0) = 2*(pow(eta(0,0),2) + pow(eta(1,0),2)) - 1;
    R(0,1) = 2*(eta(1,0)*eta(2,0) - eta(0,0)*eta(3,0));
    R(0,2) = 2*(eta(1,0)*eta(3,0) + eta(0,0)*eta(2,0));
    R(1,0) = 2*(eta(1,0)*eta(2,0) + eta(0,0)*eta(3,0));
    R(1,1) = 2*(pow(eta(0,0),2) + pow(eta(2,0),2)) - 1;
    R(1,2) = 2*(eta(2,0)*eta(3,0) - eta(0,0)*eta(1,0));
    R(2,0) = 2*(eta(1,0)*eta(3,0) - eta(0,0)*eta(2,0));
    R(2,1) = 2*(eta(2,0)*eta(3,0) + eta(0,0)*eta(1,0));
    R(2,2) = 2*(pow(eta(0,0),2) + pow(eta(3,0),2)) - 1;

    return R;

}

Eigen::Quaterniond MathUtils::eul2quat_rad(Eigen::Matrix<double, 3, 1> euler_angles_zyx)

{

    Eigen::Quaterniond q; 
    q = Eigen::AngleAxisd(euler_angles_zyx(2, 0), Eigen::Vector3d::UnitZ())
        * Eigen::AngleAxisd(euler_angles_zyx(1, 0), Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(euler_angles_zyx(0, 0), Eigen::Vector3d::UnitX());
    
    return q; 

}

Eigen::Quaterniond MathUtils::eul2quat_deg(Eigen::Matrix<double, 3, 1> euler_angles_zyx)

{

    Eigen::Quaterniond q; 
    q = Eigen::AngleAxisd(deg2rad(euler_angles_zyx(2, 0)), Eigen::Vector3d::UnitZ())
        * Eigen::AngleAxisd(deg2rad(euler_angles_zyx(1, 0)), Eigen::Vector3d::UnitY())
        * Eigen::AngleAxisd(deg2rad(euler_angles_zyx(0, 0)), Eigen::Vector3d::UnitX());
    
    return q; 

}

double MathUtils::deg2rad(double deg) 

{

    return deg * M_PI / 180;

} 

Eigen::Matrix<double, 7, 1> MathUtils::robotStates2Pose(Eigen::MatrixXd robotStates) 

{

    assertm(robotStates.rows() == 12, "robotStates must have 12 rows, 1 col."); 
    assertm(robotStates.cols() == 1, "robotStates must have 12 rows, 1 col."); 

    Eigen::Matrix<double, 7, 1> pose; 
    Eigen::Matrix3d R; 
    Eigen::MatrixXd R_tmp(9, 1);
    Eigen::Quaterniond eta; 
    R_tmp << robotStates.block<9, 1>(3, 0); 
    R_tmp.resize(3, 3);
    R = R_tmp; 
    eta = MathUtils::rot2quat(R); 
    pose << robotStates.block<3, 1>(0, 0), eta.w(), eta.vec(); 

    return pose;

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
