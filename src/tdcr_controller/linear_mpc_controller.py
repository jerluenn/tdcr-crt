import sys
sys.path.insert(0, '../classes')
sys.path.insert(0, '../utils')

import shutil
from acados_template import AcadosModel
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from casadi import *
import numpy as np
import scipy
from scipy import linalg
import time

"""
This script provides a vanilla type of LMPC controller. 
"""

class TDCR_MPC_builder: 

    def __init__(self, options, workspace = None): 

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

        self._options = options
        self._Tf = options['Tf']
        self._N = options['N']
        self._n_tendons = options['n_tendons']
        self._n_states = options['n_states']
        self._Q_weight = options['Q_weight']
        self._Q_weight_ = options['Q_weight_']
        self._R_weight = options['R_weight']
        self._Q_weight_t = options['Q_weight_t']
        self._Q_weight_t_ = options['Q_weight_t_']
        self._q_max = options['q_max']
        self._qdot_max = options['qdot_max']
        self._name = options['name']
        self._dir_name = 'c_generated_code_' + self._name

        if workspace is not None: 
            self._workspace = workspace
        else: 
            self._workspace = "~"

        self._removeOldData()
        self._create_controller()
        self._exportData()
        self._replaceData()

    def _create_controller(self):


        model_name = 'controller' + self._name

        q_dot = SX.sym('q_dot', self._n_tendons)
        q = SX.sym('q', self._n_tendons)
        x = SX.sym('x', self._n_states)
        J = SX.sym('J', self._n_states*self._n_tendons, 1)
        x_dot = SX.sym('x_dot', self._n_states)

        states = vertcat(x, q)
        states_dot = vertcat(x_dot, q_dot)
        xdot = reshape(J, self._n_states, self._n_tendons)@q_dot
        f_expl_expr = vertcat(xdot, q_dot)
        f_impl_expr = f_expl_expr - states_dot

        model = AcadosModel()
        model.name = self._name
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

        ocp.dims.N = self._N

        ocp.code_export_directory = self._dir_name + '/' + self._name
        self._code_export_dir = ocp.code_export_directory

        ocp.cost.cost_type = 'LINEAR_LS'
        ocp.cost.cost_type_e = 'LINEAR_LS'

        ocp.cost.W = scipy.linalg.block_diag(self._Q_weight, self._Q_weight_, self._R_weight)
        ocp.cost.W_e = scipy.linalg.block_diag(self._Q_weight_t, self._Q_weight_t_)

        ocp.cost.Vx = np.zeros((ny, nx))
        ocp.cost.Vx[:nx, :nx] = np.eye(nx)

        Vu = np.zeros((ny, nu))
        Vu[nx:, :] = np.eye(nu)
        ocp.cost.Vu = Vu

        ocp.cost.Vx_e = np.eye(nx)

        ocp.cost.yref = np.zeros((ny, ))
        ocp.cost.yref_e = np.zeros(ny_e, )

        ocp.constraints.idxbu = np.array([i for i in range(nu)])
        ocp.constraints.lbu = np.array([-self._qdot_max for i in range(nu)])
        ocp.constraints.ubu = np.array([self._qdot_max for i in range(nu)])
        ocp.constraints.x0 = np.zeros(nx)

        ocp.constraints.idxbx = np.array([i+self._n_states for i in range(nx - self._n_states)])
        ocp.constraints.lbx = np.array([0 for i in range(nx - self._n_states)])
        ocp.constraints.ubx = np.array([self._q_max for i in range(nx - self._n_states)])

        ocp.solver_options.levenberg_marquardt = 0.0
        # ocp.solver_options.regularize_method = 'CONVEXIFY'
        ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM'  # FULL_CONDENSING_QPOASES
        # ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'  # FULL_CONDENSING_QPOASES
        ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
        ocp.solver_options.integrator_type = 'ERK'
        ocp.solver_options.nlp_solver_type = 'SQP_RTI'  # SQP_RTI
        ocp.solver_options.nlp_solver_max_iter = 200
        ocp.parameter_values = np.zeros(self._n_states*self._n_tendons)
        # ocp.solver_options.qp_solver_cond_N = self._N

        # set prediction horizon
        ocp.solver_options.tf = self._Tf

        acados_ocp_solver = AcadosOcpSolver(
            ocp, json_file='acados_ocp_' + model.name + '.json')
        acados_integrator = AcadosSimSolver(
            ocp, json_file='acados_ocp_' + model.name + '.json')

        return acados_ocp_solver, acados_integrator


    def _replaceData(self): 

        textToReplace = "struct " + self._name + "_solver_capsule"
        textToReplace2 = "} " + self._name + "_solver_capsule"

        replacedText = textToReplace + "_replaced"
        replacedText2 = textToReplace2 + "_replaced"

        os.chdir(os.path.expanduser(
            self._workspace + "/tdcr_crt/src/tdcr_controller/" + self._dir_name + "/" + self._name))
        fullFileName = "acados_solver_" + self._name

        # Read in the file
        with open(fullFileName + '.h', 'r') as file:
            filedata = file.read()

        # Replace the target string
        filedata = filedata.replace(textToReplace, replacedText)
        filedata = filedata.replace(textToReplace2, replacedText2)

        # Write the file out again
        with open(fullFileName + '.h', 'w') as file:
            file.write(filedata)

    def _removeOldData(self):

        os.chdir(os.path.expanduser(self._workspace + "/tdcr_crt/lib/shared"))
        try: 
            os.chdir(os.path.expanduser(
            self._workspace + "/tdcr_crt/src/tdcr_controller/" + self._dir_name))
        except: 
            print("No such directory as ~/tdcr_crt/src/tdcr_controller/ + self._dir_name, changing directory to ~/tdcr_crt/src/tdcr_controller/.")
            os.chdir(os.path.expanduser(
            self._workspace + "/tdcr_crt/src/tdcr_controller/"))
        
        list_dir = [x[0] for x in os.walk(os.getcwd())]

        if os.getcwd() == os.path.expanduser(self._workspace + "/tdcr_crt/src/tdcr_controller/" + self._dir_name):

            for file_name in os.listdir(os.getcwd()):
                # construct full file path
                file = os.getcwd() + file_name
                if os.path.isfile(file):
                    print('Deleting file:', file)
                    os.remove(file)

            for _dir in list_dir:

                try:

                    print("Deleting " + _dir)
                    shutil.rmtree(_dir)

                except:

                    pass

            os.chdir("..")

        else:

            print("Removing files in the wrong folder!")

        

    def _exportData(self): 

        os.chdir(os.path.expanduser(self._workspace + "/tdcr_crt/src/tdcr_controller"))
        os.chdir(self._dir_name)
        os.chdir(self._name)
        os.system("mv *.so ../../../../lib/shared")

if __name__ == "__main__":


    # options = {} 
    # options['Tf'] = 1.0
    # options['N'] = 100
    # options['n_tendons'] = 6
    # options['n_states'] = 10
    # options['Q_weight'] = linalg.block_diag(100e3 * np.eye(2), 500e3 * np.eye(2) , 1000e3 * np.eye(2), 1e-3 * np.eye(4))
    # options['Q_weight_'] = linalg.block_diag(0.1e1 * np.eye(3), 0.1e1 * np.eye(3))
    # options['R_weight'] = 1e-2 * np.eye(6)
    # options['Q_weight_t'] = 30 * linalg.block_diag(100e3 * np.eye(2), 1000e3 * np.eye(4), 1e3 * np.eye(4))
    # options['Q_weight_t_'] = 0.5e1 * np.eye(6)
    # options['q_max'] = 40.0
    # options['qdot_max'] = 30.0
    # options['name'] = 'tdcr_lmpc'


    # ENDEFFECTOR CONTROL

    options = {} 
    options['Tf'] = 1.0
    options['N'] = 100
    options['n_tendons'] = 6
    options['n_states'] = 10
    options['Q_weight'] = linalg.block_diag(1e-5 * np.eye(2), 1e-5 * np.eye(2) , 1000e3 * np.eye(2), 100e3 * np.eye(4))
    options['Q_weight_'] = linalg.block_diag(0.1e1 * np.eye(3), 0.1e1 * np.eye(3))
    options['R_weight'] = 1e-2 * np.eye(6)
    options['Q_weight_t'] = 30 * linalg.block_diag(1e-5 * np.eye(4), 1000e3 * np.eye(2), 100e3 * np.eye(4))
    options['Q_weight_t_'] = 0.5e1 * np.eye(6)
    options['q_max'] = 40.0
    options['qdot_max'] = 30.0
    options['name'] = 'tdcr_lmpc'


    # options = {} 
    # options['Tf'] = 1.0
    # options['N'] = 100
    # options['n_tendons'] = 6
    # options['n_states'] = 6
    # options['Q_weight'] = linalg.block_diag(100e4 * np.eye(2), 100e2 * np.eye(4))
    # options['Q_weight_'] = linalg.block_diag(0.1e1 * np.eye(1), 0.3e1 * np.eye(1), 0.1e1 * np.eye(1), 0.3e1 * np.eye(1), 1e1 * np.eye(1), 0.1e1 * np.eye(1))
    # options['R_weight'] = 1e-5 * np.eye(6)
    # options['Q_weight_t'] = 30 * linalg.block_diag(100e4 * np.eye(2), 100e2 * np.eye(4))
    # options['Q_weight_t_'] = 0.5e1 * np.eye(6)
    # options['q_max'] = 40.0
    # options['qdot_max'] = 30.0
    # options['name'] = 'tdcr_lmpc'

    # TDCR_MPC_builder(options)
    TDCR_MPC_builder(options, sys.argv[1])