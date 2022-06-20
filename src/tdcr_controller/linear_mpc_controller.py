import sys
sys.path.insert(0, '../classes')
sys.path.insert(0, '../utils')

from acados_template import AcadosModel
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from casadi import *
import numpy as np
import scipy
import time

"""
This script provides a vanilla type of LMPC controller. 
"""

class MPC_builder: 

    def __init__(self, options): 

        self._options = options

    def _create_controller(Tf, N, n_tendons, n_states, Q_weight, Q_weight_, R_weight, Q_weight_t, Q_weight_t_, q_max, qdot_max, name):

        """
        Parameters
        __________

        Tf: Final time (s).
        N: Number of steps in horizon. 
        n_tendons: Number of actuation inputs, tendons. 
        n_states: Number of CARTESIAN/ORIENTATION states to be controlled. Total number of states = n_states + n_tendons
        Q_weight: n_states*n_states penalty matrix on cartersian/orientation states.
        Q_weight_: n_tendons*n_tendons penalty weight on actuation states.
        R_weight: n_tendons*n_tendons penalty weight on control inputs. 
        Q_weight_t: n_states*n_states terminal penalty matrix on cartersian/orientation states. 
        Q_weight_t_: n_states*n_states terminal penalty matrix on actuation states
        q_max: Maximum tension in each tendons. 
        qdot_max: Maximum rate of change in tendon tension.
        name: name of mpc.
        """

        model_name = 'controller' + name

        q_dot = SX.sym('q_dot', n_tendons)
        q = SX.sym('q', n_tendons)
        x = SX.sym('x', n_states)
        J = SX.sym('J', n_states*n_tendons, 1)
        x_dot = SX.sym('x_dot', n_states)

        states = vertcat(x, q)
        states_dot = vertcat(x_dot, q_dot)
        xdot = reshape(J, n_states, n_tendons)@q_dot
        f_expl_expr = vertcat(xdot, q_dot)
        f_impl_expr = f_expl_expr - states_dot

        model = AcadosModel()
        model.name = model_name
        model.x = states
        model.xdot = states_dot
        model.u = q_dot
        model.f_expl_expr = f_expl_expr
        model.f_impl_expr = f_impl_expr
        model.z = []
        model.p = J

        ocp = AcadosOcp()
        ocp.model = model

        nx = model.x.size()[0]
        nu = model.u.size()[0]
        ny = nx + nu
        ny_e = nx

        ocp.dims.N = N

        ocp.cost.cost_type = 'LINEAR_LS'
        ocp.cost.cost_type_e = 'LINEAR_LS'

        ocp.cost.W = scipy.linalg.block_diag(Q_weight, Q_weight_, R_weight)
        ocp.cost.W_e = scipy.linalg.block_diag(Q_weight_t, Q_weight_t_)

        ocp.cost.Vx = np.zeros((ny, nx))
        ocp.cost.Vx[:nx, :nx] = np.eye(nx)

        Vu = np.zeros((ny, nu))
        Vu[nx:, :] = np.eye(nu)
        ocp.cost.Vu = Vu

        ocp.cost.Vx_e = np.eye(nx)

        ocp.cost.yref = np.zeros((ny, ))
        ocp.cost.yref_e = np.zeros(ny_e, )

        ocp.constraints.idxbu = np.array([i for i in range(nu)])
        ocp.constraints.lbu = np.array([-qdot_max for i in range(nu)])
        ocp.constraints.ubu = np.array([qdot_max for i in range(nu)])
        ocp.constraints.x0 = np.zeros(nx)

        ocp.constraints.idxbx = np.array([i+n_states for i in range(nx - n_states)])
        ocp.constraints.lbx = np.array([0 for i in range(nx - n_states)])
        ocp.constraints.ubx = np.array([q_max for i in range(nx - n_states)])

        # ocp.solver_options.levenberg_marquardt = 0.001
        # ocp.solver_options.regularize_method = 'CONVEXIFY'
        ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'  # FULL_CONDENSING_QPOASES
        # ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'  # FULL_CONDENSING_QPOASES
        ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
        ocp.solver_options.integrator_type = 'ERK'
        ocp.solver_options.nlp_solver_type = 'SQP_RTI'  # SQP_RTI
        ocp.parameter_values = np.zeros(n_states*n_tendons)
        ocp.solver_options.qp_solver_cond_N = N

        # set prediction horizon
        ocp.solver_options.tf = Tf

        acados_ocp_solver = AcadosOcpSolver(
            ocp, json_file='acados_ocp_' + model.name + '.json')
        acados_integrator = AcadosSimSolver(
            ocp, json_file='acados_ocp_' + model.name + '.json')

        return acados_ocp_solver, acados_integrator


    def _removeOldData(self):

        pass 

    def _exportData(self): 

        pass 