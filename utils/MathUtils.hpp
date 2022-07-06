#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <vector>
#include <chrono>

namespace MathUtils {

    Eigen::Matrix<double, 3, 1> so3toVec(Eigen::Matrix3d so3mat);
    Eigen::Matrix<double, 6, 1> se3toVec(Eigen::Matrix4d se3mat);
    Eigen::MatrixXd forwardFiniteDifferences(Eigen::MatrixXd mat, Eigen::MatrixXd mat_plus, double eps);
    Eigen::Matrix<double, 6, 6> adjointTransformation(Eigen::Affine3d T);
    Eigen::Matrix<double, 3, 3> skew_m(Eigen::Matrix<double, 3, 1> v);
    Eigen::Matrix<double, 3, 3> skew_v(Eigen::Vector<double, 3> v);
    Eigen::Matrix3d quat2Rot(Eigen::Matrix<double, 4, 1> eta);
    Eigen::Quaterniond rot2quat(Eigen::Matrix3d rot);
    Eigen::Quaterniond eul2quat_rad(Eigen::Matrix<double, 3, 1> euler_angles_zyx);
    Eigen::Quaterniond eul2quat_deg(Eigen::Matrix<double, 3, 1> euler_angles_zyx);
    Eigen::Matrix<double, 7, 1> robotStates2Pose(Eigen::MatrixXd robotStates);
    double deg2rad(double deg);

    class Timer{

        private: 

            std::chrono::time_point<std::chrono::high_resolution_clock> t_start;
            std::chrono::time_point<std::chrono::high_resolution_clock> t_end;

        public: 

            void tic();
            void toc();
            double elapsed_time_ms;

    };

} 

#endif