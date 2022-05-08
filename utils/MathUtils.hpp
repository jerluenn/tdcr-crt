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