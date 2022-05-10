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
import rod_parameterbuilder as rpb

class TetherUnit: 

    def __init__(self, _rod_builder): 

        self._params = _rod_builder.params 
        self._dir_name = 'c_generated_code_' + 'tether_unit'
        self._buildTetherModel()

    def _buildTetherModel(self): 

        self._mass_distribution = self._params['mass_distribution']
        self._Kse = self._params['Kse']
        self._Kbt = self._params['Kbt']
        self._tether_length = self._params['tether_length']
        self._Integrator = None
        self._stepIntegrator = None
        self._integrationSteps = self._params['integration_steps']

        self._initialiseStates()
        self._createIntegrator()
        self._createstepIntegrator()

    def _initialiseStates(self):

        # Initialise all ODE states.

        self._p = SX.sym('p', 3)
        self._R = SX.sym('R', 9)
        self._n = SX.sym('n', 3)
        self._m = SX.sym('m', 3)
        self._tau = SX.sym('tau', 1)
        self._alpha = SX.sym('alpha', 1)

        self._p_d = SX.sym('p_dot', 3)
        self._R_d = SX.sym('R_dot', 9)
        self._n_d = SX.sym('n_dot', 3)
        self._m_d = SX.sym('m_dot', 3)
        self._tau_d = SX.sym('tau_dot', 1)
        self._alpha_d = SX.sym('alpha_dot', 1)
        self._Kappa_d = SX.sym('Kappa_d', 1)

        # Initialise constants

        self._g = SX([9.81, 0, 0])
        self._f_ext = self._mass_distribution * self._g
        self._Kappa = SX.sym('Kappa', 1)

        # Intermediate states

        self._u = inv(self._Kbt)@transpose(reshape(self._R, 3, 3))@self._m
        self._v = SX([0, 0, 1])


    def _createIntegrator(self):

        model_name = 'tetherunit'

        p_dot = reshape(self._R, 3, 3) @ self._v
        R_dot = reshape(self._R, 3, 3) @ skew(self._u)
        n_dot = - (self._f_ext) 
        m_dot = - cross(p_dot, self._n) 
        tau_dot = -self._Kappa*self._tau*norm_2(self._u)
        alpha_dot = 1
        kappa_dot = 0

        x = vertcat(self._p, self._R, self._n, self._m, self._tau, self._alpha, self._Kappa)
        xdot = vertcat(p_dot, reshape(R_dot, 9, 1),
                       n_dot, m_dot, tau_dot, alpha_dot, kappa_dot)
        x_dot_impl = vertcat(self._p_d, self._R_d, self._n_d, self._m_d, self._tau_d, self._alpha_d, self._Kappa_d)

        model = AcadosModel()
        model.name = model_name
        model.x = x 
        model.f_expl_expr = xdot 
        model.f_impl_expr = xdot - x_dot_impl
        model.u = SX([])
        model.z = SX([])
        model.xdot = x_dot_impl

        sim = AcadosSim()
        sim.model = model 

        Sf = self._tether_length

        sim.code_export_directory = self._dir_name + '/' + model_name
        # for exporting data to library folder afterwards.

        sim.solver_options.T = Sf
        sim.solver_options.integrator_type = 'ERK'
        sim.solver_options.num_stages = 4
        sim.solver_options.num_steps = 10

        acados_integrator = AcadosSimSolver(sim)

        self._Integrator = acados_integrator

        return acados_integrator

    def _createstepIntegrator(self):

        model_name = 'tetherunit_stepIntegrator'

        p_dot = reshape(self._R, 3, 3) @ self._v
        R_dot = reshape(self._R, 3, 3) @ skew(self._u)
        n_dot = - (self._f_ext)
        m_dot = - cross(p_dot, self._n)
        tau_dot = -self._Kappa*self._tau*norm_2(inv(self._Kbt)@(transpose(R_dot)@self._m + transpose(reshape(self._R, 3, 3))@m_dot))
        alpha_dot = 1
        kappa_dot = 0

        x = vertcat(self._p, self._R, self._n, self._m, self._tau, self._alpha, self._Kappa)
        xdot = vertcat(p_dot, reshape(R_dot, 9, 1),
                       n_dot, m_dot, tau_dot, alpha_dot, kappa_dot)

        model = AcadosModel()
        model.name = model_name
        model.x = x 
        model.f_expl_expr = xdot 
        model.u = SX([])
        model.z = []

        sim = AcadosSim()
        sim.model = model 

        Sf = self._tether_length/self._integrationSteps

        sim.code_export_directory = self._dir_name + '/' + model_name
        # for exporting data to library folder afterwards.

        sim.solver_options.T = Sf
        sim.solver_options.integrator_type = 'ERK'
        sim.solver_options.num_stages = 4
        sim.solver_options.num_steps = 1
        acados_integrator = AcadosSimSolver(sim)

        self._stepIntegrator = acados_integrator

        return acados_integrator