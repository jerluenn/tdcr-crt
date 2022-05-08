#include <IntegratorInterface.hpp>

IntegrationInterface::IntegrationInterface(sim_solver_capsule* capsule_){

    setCapsule(capsule_);

} 

IntegrationInterface::~IntegrationInterface() {};

void IntegrationInterface::setCapsule(sim_solver_capsule* capsule_arg){

    capsule = capsule_arg;

}

Eigen::MatrixXd IntegrationInterface::integrate(Eigen::MatrixXd x0) {

    sim_in_set(capsule->acados_sim_config, capsule->acados_sim_dims,
        capsule->acados_sim_in, "x", x0.data());

    status = sim_solve(capsule->acados_sim_solver, capsule->acados_sim_in, capsule->acados_sim_out);

    if (status != ACADOS_SUCCESS)
    {
        printf("acados_solve() failed with status %d.\n", status);
    }

    sim_out_get(capsule->acados_sim_config, capsule->acados_sim_dims,
        capsule->acados_sim_out, "x", x0.data());

    return x0;

}



