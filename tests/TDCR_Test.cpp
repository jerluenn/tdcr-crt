#define BOOST_LOG_DYN_LINK 1
#include <MultistageTDCR_Solver.hpp>
#include <IntegratorInterface.hpp>
#include <LevenbergMarquardtFunctor.hpp>
#include <TDCR_Interface.hpp>
#include "acados_sim_solver_multistage_straight_integrator1.h"
#include "acados_sim_solver_multistage_straight_integrator2.h"
#include "acados_sim_solver_multistage_straight_step_integrator1.h"
#include "acados_sim_solver_multistage_straight_step_integrator2.h"
#include "acados_solver_tdcr_lmpc.h"
#include "ros/ros.h"

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
    tdcr_lmpc_solver_capsule *nlpcapsule = tdcr_lmpc_acados_create_capsule();
    tdcr_lmpc_acados_create(nlpcapsule);

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
    double radius = 0.035;

    for (int i = 0; i < 3; ++i) {

        routing(0, i) = radius*cos(angle); 
        routing(1, i) = radius*sin(angle);
        routing(0, i + 3) = radius*cos(angle2);
        routing(1, i + 3) = radius*sin(angle2);
        angle += 2*PI/3;
        angle2 += 2*PI/3;

    } 

    MultistageTDCR_Solver tendon_robot(20, 6, 2, i, is, stage_tendons, routing);
    ControllerInterface controller(nlpcapsule);
    TDCR_Interface c(&tendon_robot, &controller); 

    Eigen::MatrixXd tau(6, 1); 
    tau << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; 
    c.solveForwardKinematics(tau, false);

    Eigen::MatrixXi w1(2, 1), w2(6, 1);
    Eigen::MatrixXd desiredPose(8, 1), controlInput; 
    std::vector<Eigen::MatrixXi> CSM; 
    Eigen::MatrixXi CS(2, 1); 
    CS << 0, 1;
    w1 << 0, 1; 
    w2 << 0, 1, 3, 4, 5, 6; 
    CSM.push_back(w1); 
    CSM.push_back(w2); 
    c.setDimensions(8, CSM, CS);

    desiredPose << 0.0, 0.04, -0.2, 0.2, 1., 0., 0., 0.;

    c.timer.tic();

    for (int i = 0; i < 500; ++i) 
    
    {

        controlInput = c.getHighLevelControl(desiredPose);
        
        c.simulateStep(controlInput);
        // std::cout << "Checking Boundary Conditions: " << c.checkBoundaryConditions() << "\n"; 
       

    }

    c.timer.toc();

    c.saveFullStates("test.csv");
    c.saveFullStates("test.csv");

    // std::cout << "Error: \n" << c.getCustomPoseError() << "\n\n";
    std::cout << "customPose: " << c.getCustomPose().transpose() << std::endl;
    std::cout << "Required Tau: \n" << c.getTau().transpose() << std::endl;

    free(capsule1);
    free(capsule1step);
    free(capsule2);
    free(capsule2step);

    return 0; 

}   