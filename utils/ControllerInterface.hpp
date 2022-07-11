#ifndef CONTROLLERINFERFACE_H
#define CONTROLLERINTERFACE_H

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <Eigen/Dense>
#include <MathUtils.hpp>

// acados
#include "acados/utils/print.h"
#include "acados/utils/math.h"
#include "acados_c/sim_interface.h"
#include "acados/utils/timing.h"

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"

#include "acados_c/sim_interface.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

#define assertm(exp, msg) assert(((void)msg, exp))

typedef struct tdcr_lmpc_solver_capsule
{
    // acados objects
    ocp_nlp_in *nlp_in;
    ocp_nlp_out *nlp_out;
    ocp_nlp_solver *nlp_solver;
    void *nlp_opts;
    ocp_nlp_plan *nlp_solver_plan;
    ocp_nlp_config *nlp_config;
    ocp_nlp_dims *nlp_dims;

    // number of expected runtime parameters
    unsigned int nlp_np;

    /* external functions */
    // dynamics
    external_function_param_casadi *forw_vde_casadi;
    external_function_param_casadi *expl_ode_fun;
    external_function_param_casadi *hess_vde_casadi;
    external_function_param_casadi *impl_dae_fun;
    external_function_param_casadi *impl_dae_fun_jac_x_xdot_z;
    external_function_param_casadi *impl_dae_jac_x_xdot_u_z;
    external_function_param_casadi *impl_dae_hess;
    external_function_param_casadi *gnsf_phi_fun;
    external_function_param_casadi *gnsf_phi_fun_jac_y;
    external_function_param_casadi *gnsf_phi_jac_y_uhat;
    external_function_param_casadi *gnsf_f_lo_jac_x1_x1dot_u_z;
    external_function_param_casadi *gnsf_get_matrices_fun;
    external_function_param_casadi *discr_dyn_phi_fun;
    external_function_param_casadi *discr_dyn_phi_fun_jac_ut_xt;
    external_function_param_casadi *discr_dyn_phi_fun_jac_ut_xt_hess;

    // cost
    external_function_param_casadi *cost_y_fun;
    external_function_param_casadi *cost_y_fun_jac_ut_xt;
    external_function_param_casadi *cost_y_hess;
    external_function_param_casadi *ext_cost_fun;
    external_function_param_casadi *ext_cost_fun_jac;
    external_function_param_casadi *ext_cost_fun_jac_hess;

    external_function_param_casadi cost_y_0_fun;
    external_function_param_casadi cost_y_0_fun_jac_ut_xt;
    external_function_param_casadi cost_y_0_hess;
    external_function_param_casadi ext_cost_0_fun;
    external_function_param_casadi ext_cost_0_fun_jac;
    external_function_param_casadi ext_cost_0_fun_jac_hess;

    external_function_param_casadi cost_y_e_fun;
    external_function_param_casadi cost_y_e_fun_jac_ut_xt;
    external_function_param_casadi cost_y_e_hess;
    external_function_param_casadi ext_cost_e_fun;
    external_function_param_casadi ext_cost_e_fun_jac;
    external_function_param_casadi ext_cost_e_fun_jac_hess;

    // constraints
    external_function_param_casadi *phi_constraint;
    external_function_param_casadi *nl_constr_h_fun_jac;
    external_function_param_casadi *nl_constr_h_fun;
    external_function_param_casadi *nl_constr_h_fun_jac_hess;

    external_function_param_casadi phi_e_constraint;
    external_function_param_casadi nl_constr_h_e_fun_jac;
    external_function_param_casadi nl_constr_h_e_fun;
    external_function_param_casadi nl_constr_h_e_fun_jac_hess;
} nlp_solver_capsule;

class ControllerInterface 

{

    public: 

        ControllerInterface(tdcr_lmpc_solver_capsule* capsule_); 
        virtual ~ControllerInterface();
        Eigen::MatrixXd solveOptimalControl(Eigen::MatrixXd currentPose, Eigen::MatrixXd J, Eigen::MatrixXd y); 
        MathUtils::Timer timer;

    private: 

        int status; 
        Eigen::MatrixXd J;
        Eigen::MatrixXd y;
        tdcr_lmpc_solver_capsule* capsule;
        void setCapsule(tdcr_lmpc_solver_capsule* capsule_arg);
        int setJacobians(Eigen::MatrixXd J);
        void setPose(Eigen::MatrixXd currentPose);
        void setReference(Eigen::MatrixXd y);
        int N; // Horizon
        int NU; // Num Control States
        int NX_noInput;

} ;

#endif