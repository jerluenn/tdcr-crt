#define BOOST_LOG_DYN_LINK 1
#include <MultistageTDCR_Solver.hpp>
#include <IntegratorInterface.hpp>
#include "acados_sim_solver_multistage_straight_integrator1.h"
#include "acados_sim_solver_multistage_straight_integrator2.h"
#include "acados_sim_solver_multistage_straight_step_integrator1.h"
#include "acados_sim_solver_multistage_straight_step_integrator2.h"


#define PI 3.14159 

using namespace MathUtils;

int main() {

    sim_solver_capsule *capsule1 = multistage_straight_integrator1_acados_sim_solver_create_capsule();
    multistage_straight_integrator1_acados_sim_create(capsule1);
    sim_solver_capsule *capsule2 = multistage_straight_integrator2_acados_sim_solver_create_capsule();
    multistage_straight_integrator2_acados_sim_create(capsule2);
    sim_solver_capsule *capsule1step = multistage_straight_integrator1_acados_sim_solver_create_capsule();
    multistage_straight_step_integrator1_acados_sim_create(capsule1step);
    sim_solver_capsule *capsule2step = multistage_straight_integrator2_acados_sim_solver_create_capsule();
    multistage_straight_step_integrator2_acados_sim_create(capsule2step);

    std::vector<IntegrationInterface> i;
    std::vector<IntegrationInterface> is;

    IntegrationInterface i1(capsule1), i2(capsule2), ii1(capsule1step), ii2(capsule2step); 
    i.push_back(i1);
    i.push_back(i2);
    is.push_back(ii1);
    is.push_back(ii2);

    Eigen::MatrixXd stage_tendons;
    stage_tendons.resize(2, 6);
    stage_tendons << 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1;

    Eigen::MatrixXd routing; 
    routing.resize(3, 6);
    routing.row(2).setZero(); 
    double angle = 0.0; 
    double angle2 = PI;
    double radius = 0.01506;

    for (int i = 0; i < 3; ++i) {

        routing(0, i) = radius*cos(angle); 
        routing(1, i) = radius*sin(angle);
        routing(0, i + 3) = radius*cos(angle2);
        routing(1, i + 3) = radius*sin(angle2);
        angle += 2*PI/3;
        angle2 += 2*PI/3;

    } 

    
    std::cout << routing << "\n\n\n"; 

    MultistageTDCR_Solver b(6, 2, i, is, stage_tendons, routing);

    // b.getRobotStates(true, true);

    free(capsule1);
    free(capsule1step);
    free(capsule2);
    free(capsule2step);

    Eigen::Affine3d T; 
    Eigen::Matrix3d R; 
    Eigen::Vector3d V(0, 0, 3); 
    R << 1, 0, 0, 0, 0, -1, 0, 1, 0; 
    T.linear() = R; 
    T.translation() = V; 

    std::cout << MathUtils::adjointTransformation(T) << "\n\n";


    return 0; 

}   