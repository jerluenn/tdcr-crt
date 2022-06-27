#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <vector>
#include <chrono>

#define assertm(exp, msg) assert(((void)msg, exp))

typedef Eigen::MatrixXd(*func)(Eigen::MatrixXd);

namespace Integrator {

typedef Eigen::MatrixXd(*func)(Eigen::MatrixXd);
template<func f>
Eigen::MatrixXd integrate_step_rk4_autonomous(Eigen::MatrixXd x0, double dt)

{

    assertm(x0.cols() == 1, "x0 must only have one column.");
    Eigen::MatrixXd k1(x0.rows(), 1), k2(x0.rows(), 1), k3(x0.rows(), 1), k4(x0.rows(), 1), output(x0.rows(), 1);
    k1 = f(x0); 
    k2 = f(x0 + (dt/2)*k1);
    k3 = f(x0 + (dt/2)*k2);
    k4 = f(x0 + dt*k3);

    output = x0 + (1.0/6.0) * dt * (k1 + 2*k2 + 2*k3 + k4); 

    return output; 

}


}

#endif
