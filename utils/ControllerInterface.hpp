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
    ocp_nlp_out *sens_out;
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
    

} tdcr_lmpc_solver_capsule;

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