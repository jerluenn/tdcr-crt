#include <ControllerInterface.hpp>

ControllerInterface::ControllerInterface(nlp_solver_capsule* capsule_) 

{

    setCapsule(capsule_);
    
}

ControllerInterface::~ControllerInterface() {};

void ControllerInterface::setCapsule(nlp_solver_capsule* capsule_arg) 

{

    capsule = capsule_arg;
    N = capsule->nlp_dims->N;
    NU = *capsule->nlp_dims->nu;
    NX_noInput = *capsule->nlp_dims->nx - NU;

}

Eigen::MatrixXd ControllerInterface::solveOptimalControl(Eigen::MatrixXd currentPose, Eigen::MatrixXd J,  Eigen::MatrixXd y) 

{

    Eigen::MatrixXd optimalControl; 
    optimalControl.resize(NU, 1);
    Eigen::MatrixXd states(14, 1);
    states.setZero();

    setJacobians(J);
    setPose(currentPose); 
    setReference(y);

    status = ocp_nlp_solve(capsule->nlp_solver, capsule->nlp_in, capsule->nlp_out);

    if (status != ACADOS_SUCCESS)
    {
        printf("acados_solve() failed with status %d.\n", status);
    }

    ocp_nlp_out_get(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_out, 0, "u", optimalControl.data());
    ocp_nlp_out_get(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_out, 0, "x", states.data());

    return optimalControl;

}

void ControllerInterface::setPose(Eigen::MatrixXd currentPose) 

{

    status = ocp_nlp_constraints_model_set(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_in, 0, "lbx", currentPose.data());
    status = ocp_nlp_constraints_model_set(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_in, 0, "ubx", currentPose.data()); 

}

void ControllerInterface::setReference(Eigen::MatrixXd y) 

{

    assertm(y.cols() == 1, "y must be a vector.");
    

    for (int i = 0; i < N; ++i)

    {

        status = ocp_nlp_cost_model_set(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_in, i, "yref", y.data());

    }


}

int ControllerInterface::setJacobians(Eigen::MatrixXd J)
{
    int solver_status = 0;
    int np = J.rows()*J.cols();

    int casadi_np = NU*NX_noInput;
    if (casadi_np != np) {
        printf("acados_update_params: trying to set %i parameters for external functions."
            " External function has %i parameters. Exiting.\n", np, casadi_np);
        exit(1);
    }

    for (int stage = 0; stage < N; ++stage) 
    
    {

        if (stage < N && stage >= 0)
        {
            capsule->forw_vde_casadi[stage].set_param(capsule->forw_vde_casadi+stage, J.data());
            capsule->expl_ode_fun[stage].set_param(capsule->expl_ode_fun+stage, J.data());
        

            // constraints

            // cost
            if (stage == 0)
            {
            }
            else // 0 < stage < N
            {
            }
        }

        else // stage == N
        {
            // terminal shooting node has no dynamics
            // cost
            // constraints
        
        }

    }

    return solver_status;
}

