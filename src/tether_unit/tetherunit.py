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
        self._eta = SX.sym('self._eta', 4) 
        # self._R = SX.sym('R', 9)
        self._n = SX.sym('n', 3)
        self._m = SX.sym('m', 3)
        self._tau = SX.sym('tau', 1)
        self._alpha = SX.sym('alpha', 1)

        self._p_d = SX.sym('p_dot', 3)
        self._eta_d = SX.sym('eta_dot', 4)
        self._n_d = SX.sym('n_dot', 3)
        self._m_d = SX.sym('m_dot', 3)
        self._tau_d = SX.sym('tau_dot', 1)
        self._alpha_d = SX.sym('alpha_dot', 1)
        self._Kappa_d = SX.sym('Kappa_d', 1)

        # Initialise constants

        self._g = SX([9.81, 0, 0])
        self._f_ext = self._mass_distribution * self._g
        self._Kappa = SX.sym('Kappa', 1)

        # Setting R 

        self._R = SX(3,3)
        self._R[0,0] = 2*(self._eta[0]**2 + self._eta[1]**2) - 1
        self._R[0,1] = 2*(self._eta[1]*self._eta[2] - self._eta[0]*self._eta[3])
        self._R[0,2] = 2*(self._eta[1]*self._eta[3] + self._eta[0]*self._eta[2])
        self._R[1,0] = 2*(self._eta[1]*self._eta[2] + self._eta[0]*self._eta[3])
        self._R[1,1] = 2*(self._eta[0]**2 + self._eta[2]**2) - 1
        self._R[1,2] = 2*(self._eta[2]*self._eta[3] - self._eta[0]*self._eta[1])
        self._R[2,0] = 2*(self._eta[1]*self._eta[3] - self._eta[0]*self._eta[2])
        self._R[2,1] = 2*(self._eta[2]*self._eta[3] + self._eta[0]*self._eta[1])
        self._R[2,2] = 2*(self._eta[0]**2 + self._eta[3]**2) - 1

        # Intermediate states

        self._u = inv(self._Kbt)@transpose(reshape(self._R, 3, 3))@self._m
        self._v = SX([0, 0, 1])
        self._k = 0.5


    def _createIntegrator(self):

        model_name = 'tetherunit'

        c = self._k*(1-transpose(self._eta)@self._eta)

        p_dot = reshape(self._R, 3, 3) @ self._v
        eta_dot = vertcat(
            0.5*(-self._u[0]*self._eta[1] - self._u[1]*self._eta[2] - self._u[2]*self._eta[3]),
            0.5*(self._u[0]*self._eta[0] + self._u[2]*self._eta[2] - self._u[1]*self._eta[3]),
            0.5*(self._u[1]*self._eta[0] - self._u[2]*self._eta[1] + self._u[0]*self._eta[3]),
            0.5*(self._u[2]*self._eta[0] + self._u[1]*self._eta[1] - self._u[0]*self._eta[2])
        ) + c * self._eta 
        n_dot = - (self._f_ext) 
        m_dot = - cross(p_dot, self._n) 
        tau_dot = -self._Kappa*self._tau*norm_2(self._u)
        alpha_dot = 1
        kappa_dot = 0

        x = vertcat(self._p, self._eta, self._n, self._m, self._tau, self._alpha, self._Kappa)
        xdot = vertcat(p_dot, eta_dot,
                       n_dot, m_dot, tau_dot, alpha_dot, kappa_dot)
        x_dot_impl = vertcat(self._p_d, self._eta_d, self._n_d, self._m_d, self._tau_d, self._alpha_d, self._Kappa_d)

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
        sim.solver_options.num_steps = 20

        acados_integrator = AcadosSimSolver(sim)

        self._Integrator = acados_integrator
        self._linearisedEquations = jacobian(xdot, x)
        self.linearisedEquationsFunction = Function('linear_eq', [x], [self._linearisedEquations])

        return acados_integrator

    def _createstepIntegrator(self):

        model_name = 'tetherunit_stepIntegrator'

        k = 0.01
        c = self._k*(1-transpose(self._eta)@self._eta)

        p_dot = reshape(self._R, 3, 3) @ self._v
        eta_dot = vertcat(
            0.5*(-self._u[0]*self._eta[1] - self._u[1]*self._eta[2] - self._u[2]*self._eta[3]),
            0.5*(self._u[0]*self._eta[0] + self._u[2]*self._eta[2] - self._u[1]*self._eta[3]),
            0.5*(self._u[1]*self._eta[0] - self._u[2]*self._eta[1] + self._u[0]*self._eta[3]),
            0.5*(self._u[2]*self._eta[0] + self._u[1]*self._eta[1] - self._u[0]*self._eta[2])
        ) + c * self._eta 
        n_dot = - (self._f_ext)
        m_dot = - cross(p_dot, self._n) 
        tau_dot = -self._Kappa*self._tau*norm_2(self._u)
        alpha_dot = 1
        kappa_dot = 0

        x = vertcat(self._p, self._eta, self._n, self._m, self._tau, self._alpha, self._Kappa)
        xdot = vertcat(p_dot, eta_dot,
                       n_dot, m_dot, tau_dot, alpha_dot, kappa_dot)
        x_dot_impl = vertcat(self._p_d, self._eta_d, self._n_d, self._m_d, self._tau_d, self._alpha_d, self._Kappa_d)

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

        Sf = self._tether_length/self._integrationSteps

        sim.code_export_directory = self._dir_name + '/' + model_name
        # for exporting data to library folder afterwards.

        sim.solver_options.T = Sf
        sim.solver_options.integrator_type = 'ERK'
        sim.solver_options.num_stages = 4
        sim.solver_options.num_steps = 10
        acados_integrator = AcadosSimSolver(sim)

        self._stepIntegrator = acados_integrator

        return acados_integrator