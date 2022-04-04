#include <MultistageTDCR_Solver.hpp>
#include <IntegratorInterface.hpp>
#include "acados_sim_solver_multistage_straight_integrator1.h"
#include "acados_sim_solver_multistage_straight_integrator2.h"
#include "acados_sim_solver_multistage_straight_step_integrator1.h"
#include "acados_sim_solver_multistage_straight_step_integrator2.h"



int main() {

    sim_solver_capsule *capsule1 = multistage_straight_integrator1_acados_sim_solver_create_capsule();
    multistage_straight_integrator1_acados_sim_create(capsule1);
    sim_solver_capsule *capsule2 = multistage_straight_integrator1_acados_sim_solver_create_capsule();
    multistage_straight_integrator1_acados_sim_create(capsule2);
    sim_solver_capsule *capsule1step = multistage_straight_integrator1_acados_sim_solver_create_capsule();
    multistage_straight_integrator1_acados_sim_create(capsule1step);
    sim_solver_capsule *capsule2step = multistage_straight_integrator1_acados_sim_solver_create_capsule();
    multistage_straight_integrator1_acados_sim_create(capsule2step);

    std::vector<IntegrationInterface> i;
    std::vector<IntegrationInterface> is;

    IntegrationInterface i1(capsule1), i2(capsule2), ii1(capsule1step), ii2(capsule2step); 
    i.push_back(i1);
    i.push_back(i2);
    is.push_back(ii1);
    is.push_back(ii2);

    MultistageTDCR_Solver b(4, 2, i, is);

    b.getRobotStates(true, true);
    Eigen::MatrixXd tau;
    tau.resize(4, 1); 
    tau << 5, 0, 0, 0;
    b.setTau(tau);
    b.getRobotStates(true, true);

    Eigen::MatrixXd stage_tendons;
    stage_tendons.resize(3, 4);
    stage_tendons << 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0;
    std::cout << stage_tendons.row(0).sum() << "\n";



    return 0; 

}   