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
import tetherunit, rod_parameterbuilder
from matplotlib import pyplot as plt

class TetherUnitBoundarySolver: 

    def __init__(self, robot_dict, initConditions, distalPose): 

        self.build_tetherObject(robot_dict)
        self.boundary_length = robot_dict['tether_length'] # used for plotting
        self.integration_steps = robot_dict['integration_steps']
        self.distalPose = distalPose 
        self.initConditions = initConditions
        self.distalConditions = None 

    def build_tetherObject(self, robot_dict): 

        builder = rod_parameterbuilder.Rod_Parameter_Builder()
        builder.createHollowRod(robot_dict)
        self.tetherObject = tetherunit.TetherUnit(builder)

    def set_and_integrate(self): 

        self.tetherObject._Integrator.set('x', self.initConditions)
        self.tetherObject._Integrator.solve()
        self.distalConditions = self.tetherObject._Integrator.get('x')

        return self.distalConditions

    def solveBVP(self, debug, plot): 

        t_start = time.time()
        sol = least_squares(self.getResiduals, np.zeros(6), method = 'trf')
        self.initConditions[12:18] = sol.x 
        self.set_and_integrate()

        if debug == True: 

            print(f"Internal forces and moments at proximal end: {sol.x}.")
            print(f"Internal forces and moments at distal end: {self.distalConditions[12:18]}")
            print(f"Pose at distal end: {self.distalConditions[0: 3]}")
            print(f"Tension at proximal end: {self.initConditions[18]}")
            print(f"Tension at distal end: {self.distalConditions[18]}")
            print(f"Total cost: {sol.cost}")
            print(f"Time taken: {time.time() - t_start} s.")
            

        if plot == True: 

            self.plotData(debug)

    def plotData(self, debug):
 
        self.poseData = np.zeros((self.integration_steps + 1, 3))
        self.internalWrenchData = np.zeros((self.integration_steps + 1, 6))
        self.arcData = np.zeros((self.integration_steps + 1, 1))
        self.poseData[0, :] = self.initConditions[0:3]
        self.internalWrenchData[0, :] = self.initConditions[12:18]
        states_i = self.initConditions

        for i in range(self.integration_steps): 
            
            self.tetherObject._stepIntegrator.set('x', states_i)
            self.tetherObject._stepIntegrator.solve()
            states_i = self.tetherObject._stepIntegrator.get('x')
            self.poseData[i + 1, :] = states_i[0:3]
            self.arcData[i + 1, :] = states_i[19]
            self.internalWrenchData[i + 1, :] = states_i[12:18]

            if debug == True: 

                pass
                # print("States at node", i, ": ", states_i)

        if debug == True: 

            for i in range(6): 

                plt.plot(self.arcData, self.internalWrenchData[:,i])
                plt.show()

            

        ax = plt.axes(projection='3d')
        ax.set_xlim3d([0, self.boundary_length])
        ax.set_xlabel('Z')

        ax.set_ylim3d([-self.boundary_length/4, self.boundary_length/4])
        ax.set_ylabel('Y')

        ax.set_zlim3d([0, self.boundary_length])
        ax.set_zlabel('X')
        ax.plot3D(self.poseData[:, 2], self.poseData[:, 1], -self.poseData[:, 0])
        plt.show()

        

    def getResiduals(self, guess): 

        self.initConditions[12:18] = guess 
        self.set_and_integrate()
        residual = np.hstack((50*self.distalConditions[0:12] - 50*self.distalPose, self.initConditions[12:18]))
        # print(norm_2(residual))

        return residual

def test_Function(testClass, initCondtions): 

    testClass.initConditions = initCondtions
    print(testClass.set_and_integrate())
    print(testClass.tetherObject._Kbt)

def test_Function2(testClass, initConditions):

    testClass.initConditions = initConditions

    for i in range(5): 

        testClass.set_and_integrate()
        testClass.plotData(True)
        plt.show()
        testClass.initConditions = testClass.distalConditions
        

if __name__ == "__main__":

    robot_dict = {}
    robot_dict['type'] = 'hollow_rod'
    robot_dict['outer_radius'] = 0.002
    robot_dict['inner_radius'] = 0.0006
    robot_dict['elastic_modulus'] = 2e9
    robot_dict['mass_distribution'] = 0.035
    robot_dict['tether_length'] = 1.5
    robot_dict['shear_modulus'] = 70e9
    robot_dict['integration_steps'] = 50

    initConditions = np.array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 5, 0, 0.05])
    distalPose = np.array([-0.8, 0, 1.2, 1, 0, 0, 0, 1, 0, 0, 0, 1])
    testClass = TetherUnitBoundarySolver(robot_dict, initConditions, distalPose)
    testClass.solveBVP(True, True)

    initConditions = np.array([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 5, 0, 0.05])

    # test_Function(testClass, initConditions)
    # test_Function2(testClass, testClass.initConditions)


