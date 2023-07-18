import sys
sys.path.insert(0, '../../utils')
import os
import shutil

from acados_template import AcadosModel
from acados_template import AcadosOcp, AcadosOcpSolver, AcadosSimSolver, AcadosSim
from casadi import *

import numpy as np
from scipy.integrate import odeint
from scipy.optimize._lsq.least_squares import least_squares

import time

from tdcr_baseclass import Tendon_Robot_Builder, Robot_Parameter_Builder, Robot_Type_Builder
import bokeh


class MultistageTDCR(Tendon_Robot_Builder):

    def __init__(self, _type_params, _phyiscal_params, workspace = None):

        self._type_params = _type_params.params
        self._phyiscal_params = _phyiscal_params.params
        if workspace is not None: 
            self._workspace = workspace
        else: 
            self._workspace = "~"
        self._buildRobotModel()

    def _buildRobotModel(self):

        print("======================================================================================================================================================================================")
        print("========================================================================= Building Robot... ==========================================================================================")
        print("======================================================================================================================================================================================")

        self._r = self._type_params['routing']
        self._robot_type = self._type_params['robot_type']
        self._mass_distribution = self._phyiscal_params['mass_distribution']
        self._tip_weight = self._phyiscal_params['tip_weight']
        self._Kse = self._phyiscal_params['Kse']
        self._Kbt = self._phyiscal_params['Kbt']
        self._robot_length = self._phyiscal_params['robot_length']
        self._num_tendons = self._type_params['num_tendons']
        # vector of lengths.
        self._stage_lengths = self._type_params['stage_lengths']
        # indices of cables active in every stage.
        self._stage_tendons = self._type_params['stage_tendons']
        # total number of stages.
        self._num_stage = np.size(self._stage_lengths)
        self._N = self._type_params['integration_steps']

        ### Attritubes for export. ###

        self._dir_list = []
        self._integrator_names = []

        try:

            self._removeOldData()

        except:

            pass

        self._initialiseStates()
        self._createIntegrators()
        self._exportData()
        self._replaceData()

        print("=====================================================================================================================================================================================")
        print("====================================================================================== Build Success. ===============================================================================")
        print("=====================================================================================================================================================================================")

    def _removeOldData(self):

        os.chdir(os.path.expanduser(self._workspace + "/tdcr_crt/lib/shared"))

        os.system("rm libacados_sim_solver_multistage_straight*.so")

        self._dir_name = 'c_generated_code_' + self._robot_type
        os.chdir(os.path.expanduser(
            self._workspace + "/tdcr_crt/src/tdcr_model/" + self._dir_name))

        list_dir = [x[0] for x in os.walk(os.getcwd())]
        list_dir.pop(0)

        if os.getcwd() == os.path.expanduser(self._workspace + "/tdcr_crt/src/tdcr_model/" + self._dir_name):

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

            raise ValueError("Removing files in the wrong folder!")

        

    def _initialiseStates(self):

        # Initialise all ODE states.

        self._p = SX.sym('p', 3)
        self._R = SX.sym('R', 9)
        self._n = SX.sym('n', 3)
        self._m = SX.sym('m', 3)
        self._tau = SX.sym('tau', self._num_tendons)
        self._alpha = SX.sym('alpha', 1)

        self._p_d = SX.sym('p_dot', 3)
        self._R_d = SX.sym('R_dot', 9)
        self._n_d = SX.sym('n_dot', 3)
        self._m_d = SX.sym('m_dot', 3)
        
        self._tau_d = SX.sym('tau_dot', self._num_tendons)
        self._alpha_d = SX.sym('alpha_dot', 1)

        # Initialise constants

        self._g = SX([9.81, 0, 0])
        self._f_ext = self._mass_distribution * self._g

        # Intermediate states

        self._u = inv(self._Kbt)@transpose(reshape(self._R, 3, 3))@self._m
        self._v = SX([0, 0, 1])

        # Boundary forces, moments

        self._F_t = [SX.zeros(3) for i in range(self._num_stage)]
        self._L_t = [SX.zeros(3) for i in range(self._num_stage)]

        # Distributed forces, moments

        self._f_t = [SX.zeros(3) for i in range(self._num_stage)]
        self._l_t = [SX.zeros(3) for i in range(self._num_stage)]

        # Update forces and moments for both distributed and boundary.

        self._calc_forces_moments()

    def _createModel(self, stage):

        p_dot = reshape(self._R, 3, 3) @ self._v
        R_dot = reshape(self._R, 3, 3) @ skew(self._u)
        n_dot = - (self._f_ext + self._f_t[stage])
        m_dot = - cross(p_dot, self._n)
        tau_dot = SX.zeros(self._num_tendons)
        alpha_dot = 1

        x = vertcat(self._p, self._R, self._n, self._m, self._tau, self._alpha)
        xdot = vertcat(p_dot, reshape(R_dot, 9, 1),
                       n_dot, m_dot, tau_dot, alpha_dot)

        assert x.size()[0] == xdot.size()[
            0], f"Both x and xdot must be the same size. x.size = {x.size()[0]}, xdot.size = {xdot.size()[0]}"

        return x, xdot

    def _createIntegrators(self):

        integrators_list = []
        params = {}
        params['name'] = self._robot_type + '_integrator'
        params['num_stages'] = 4
        params['num_steps'] = 12
        params['step_bool'] = False
        params_step = {}
        params_step['name'] = self._robot_type + '_step_integrator'
        params_step['num_stages'] = 4
        params_step['num_steps'] = 2
        params_step['step_bool'] = True
        # params_step['num_steps'] = 20

        for i in range(self._num_stage):

            x, xdot = self._createModel(i)
            integrators_list.append(self._createIntegrator(i, x, xdot, params))
            integrators_list.append(
                self._createIntegrator(i, x, xdot, params_step))

        return integrators_list

    def _createIntegrator(self, stage, x, xdot, params):

        model_name = params['name'] + str(stage + 1)

        model = AcadosModel()
        model.name = model_name
        model.x = x
        model.u = SX([])
        model.f_expl_expr = xdot
        model.z = []

        sim = AcadosSim()

        sim.model = model
        Sf = self._stage_lengths[stage]

        sim.code_export_directory = self._dir_name + '/' + model_name

        # for exporting data to library folder afterwards.
        os.chdir(os.path.expanduser(
            self._workspace + "/tdcr_crt/src/tdcr_model/"))
        self._dir_list.append(sim.code_export_directory)
        self._integrator_names.append(model_name)

        if params['step_bool'] == True:

            sim.solver_options.T = Sf/self._N
        
        else: 

            sim.solver_options.T = Sf
        
        sim.solver_options.integrator_type = 'ERK'
        sim.solver_options.num_stages = params['num_stages']
        sim.solver_options.num_steps = params['num_steps']
        sim.solver_options.sens_forw = False 

        acados_integrator = AcadosSimSolver(sim)

        return acados_integrator

    def _replaceData(self):

        textToReplace = "struct sim_solver_capsule"
        textToReplace2 = "} sim_solver_capsule"

        for i in range(0, len(self._dir_list)):

            replacedText = textToReplace + self._integrator_names[i]
            replacedText2 = textToReplace2 + self._integrator_names[i]

            os.chdir(os.path.expanduser(
                self._workspace + "/tdcr_crt/src/tdcr_model/" + self._dir_list[i]))
            fullFileName = "acados_sim_solver_" + self._integrator_names[i]

            # Read in the file
            # with open(fullFileName + '.c', 'r') as file :
            #    filedata = file.read()

            # Replace the target string
            # filedata = filedata.replace(textToReplace, replacedText)

            # Write the file out again
            # with open(fullFileName + '.c', 'w') as file:
            #    file.write(filedata)

            # Read in the file
            with open(fullFileName + '.h', 'r') as file:
                filedata = file.read()

            # Replace the target string
            filedata = filedata.replace(textToReplace, replacedText)
            filedata = filedata.replace(textToReplace2, replacedText2)

            # Write the file out again
            with open(fullFileName + '.h', 'w') as file:
                file.write(filedata)

    def _calc_forces_moments(self):

        num_stages = self._num_stage
        assert(np.sum(self._stage_tendons) == self._num_tendons)

        # Find out which tendons affect which stage.

        for i in range(num_stages):

            tmp_vector = np.zeros(self._num_tendons)

            for p in range(i, num_stages, 1):

                tmp_vector += self._stage_tendons[p]

            num_tendons_list = np.where(tmp_vector == 1)

            self._f_t[i] = self._calc_distributed_force(num_tendons_list)
            self._F_t[i] = self._calc_Boundary_Force(num_tendons_list)
            self._L_t[i] = self._calc_Boundary_Moment(num_tendons_list)

    def _calc_Boundary_Force(self, tendons_index):

        F_t = SX.zeros(3)

        for i in tendons_index[0]:

            F_t += - self._tau[i] @ reshape(self._R, 3, 3) @ self._v

        return F_t

    def _calc_Boundary_Moment(self, tendons_index):

        L_t = SX.zeros(3)

        for i in tendons_index[0]:

            L_t += - self._tau[i] @ cross(reshape(self._R, 3, 3)
                                          @ self._r[:, i], reshape(self._R, 3, 3) @ self._v)

        return L_t

    def _calc_distributed_force(self, tendons_index):

        f_t = SX.zeros(3)

        for i in tendons_index[0]:

            f_t += - (self._tau[i]) * ((skew(reshape(self._R, 3, 3)@self._v)) @ (skew(
                reshape(self._R, 3, 3)@self._v)))@(reshape(self._R, 3, 3)@skew(self._u)@self._v)

        return f_t

    def getInfo(self):

        print("==========================================")
        print("============== Robot Info ================")
        print(self._type_params)
        print(self._phyiscal_params)
        print("==========================================")

        return self._type_params + self._phyiscal_params

    def _exportData(self):

        os.chdir(os.path.expanduser(self._workspace + "/tdcr_crt/src/tdcr_model"))
        os.chdir(self._dir_list[0])

        for _dir in self._dir_list:

            os.chdir("../..")
            os.chdir(_dir)
            os.system("mv *.so ../../../../lib/shared")


if __name__ == "__main__":

    ########################################
    ########## TWO STAGE ROBOT #############
    ########################################

    # robot_dict = {}
    # # robot_dict['type'] = 'hollow_rod'
    # robot_dict['type'] = 'Solid Rod'
    # # robot_dict['outer_radius'] = 0.001
    # # robot_dict['inner_radius'] = 0.0008
    # robot_dict['radius'] = 0.0015
    # robot_dict['shear_modulus'] = 200e9
    # robot_dict['elastic_modulus'] = 200e9
    # robot_dict['mass_distribution'] = 0.1
    # robot_dict['num_tendons'] = 6
    
    # robot_dict['time_step'] = 0.01
    # robot_dict['robot_length'] = 0.4
    # robot_dict['tip_weight'] = 0.0

    # physical_builder = Robot_Parameter_Builder()
    # # physical_builder.createHollowRod(robot_dict)
    # physical_builder.createRod(robot_dict)

    # robot_dict = {}
    # robot_dict['num_tendons'] = 6
    # robot_dict['routing'] = transpose(SX([[0.0350, 0.0, 0], [-0.0175, 0.0303, 0], [-0.0175, -
    #                                                                                   0.0303, 0], [-0.035, -0.0, 0], [0.0175, -0.0303, 0], [0.0175, 0.0303, 0]]))

    # robot_dict['stage_lengths'] = np.array([0.5, 0.5])
    # robot_dict['stage_tendons'] = np.array(
    #     [[1, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 1]])
    # robot_dict['integration_steps'] = 20

    # type_builder = Robot_Type_Builder()
    # type_builder.createMultiStageStraight(robot_dict)

    # # c = MultistageTDCR(type_builder, physical_builder)
    # c = MultistageTDCR(type_builder, physical_builder, sys.argv[1])

    ##########################################
    ########## THREE STAGE ROBOT #############
    ##########################################

    robot_dict = {}
    robot_dict['type'] = 'hollow_rod'
    # robot_dict['type'] = 'Solid Rod'
    robot_dict['outer_radius'] = 0.001
    robot_dict['inner_radius'] = 0.0005
    # robot_dict['radius'] = 0.0015
    robot_dict['shear_modulus'] = 200e9
    robot_dict['elastic_modulus'] = 70e9
    robot_dict['mass_distribution'] = 0.278
    robot_dict['num_tendons'] = 6
    
    robot_dict['time_step'] = 0.01
    robot_dict['robot_length'] = 0.4
    robot_dict['tip_weight'] = 0.0

    physical_builder = Robot_Parameter_Builder()
    physical_builder.createHollowRod(robot_dict)
    # physical_builder.createRod(robot_dict)

    robot_dict = {}
    robot_dict['num_tendons'] = 6
    robot_dict['routing'] = transpose(SX([[0.0189, 0.012164, 0], [-0.003202, 0.0223, 0], [-0.0189, 0.012164, 0], [0.0189, -0.012164, 0], [-0.003202, -0.0223, 0], [-0.0189, -0.012164, 0]]))

    stage_length = [0.215, 0.13, 0.125]

    robot_dict['stage_lengths'] = np.array([stage_length[0], stage_length[1], stage_length[2]])
    robot_dict['stage_tendons'] = np.array(
        [[0, 0, 1, 0, 0, 1], [0, 1, 0, 0, 1, 0], [1, 0, 0, 1, 0, 0]])
    robot_dict['integration_steps'] = 20

    type_builder = Robot_Type_Builder()
    type_builder.createMultiStageStraight(robot_dict)

    # c = MultistageTDCR(type_builder, physical_builder)
    c = MultistageTDCR(type_builder, physical_builder, sys.argv[1])
