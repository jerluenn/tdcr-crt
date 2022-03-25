from abc import abstractmethod, ABC
import numpy as np
from casadi import *


class Tendon_Robot_Builder(ABC):

    @abstractmethod
    def _initialiseStates(self):

        pass

    @abstractmethod
    def _buildRobotModel(self):

        pass

    @abstractmethod
    def getInfo(self):

        pass

    @abstractmethod
    def _createIntegrator(self):

        pass

    @abstractmethod

    def exportData(self):

        pass


class Robot_Type_Builder:

    def __init__(self): 

        self._valid_build = False

    def createMultiStageStraight(self, params):

        print("Building MultiStage straight Continuum Robot...")

        if not self._valid_build:

            self.params = params
            self.params['robot_type'] = 'multistage_straight'
            self._checkDict(self.params)

        else: 

            raise ValueError("Parameters have been built before.")

        return self.params

    def createSingleStageStraight(self, params):

        print("Building SingleStage straight Continuum Robot...")

        if not self._valid_build:

            self.params = params
            self.params['robot_type'] = 'singlestage_straight'
            self._checkDict(self.params)

        else: 

            raise ValueError("Parameters have been built before.")

        return self.params

    def createHybrid(self, params):

        print("Building Hybrid Continuum Robot...")

        if not self._valid_build:

            self.params = params
            self.params['robot_type'] = 'hybrid'
            self._checkDict(self.params)

        else: 

            raise ValueError("Parameters have been built before.")

        return self.params

    def _checkDict(self, params):


        try:

            assert type(params['robot_type']) == str, "robot_type must be a str type."
            if type(params['routing']) == SX:
                self.params['num_tendons'] = self.params['routing'].size()[1]
            else:
                raise ValueError("Routing type must be 'casadi.casadi.SX'. ")
            assert type(params['num_tendons']) == int, "num_tendons must be an integer!"
            if params['robot_type'] == 'multistage_straight': 
                assert type(params['stage_lengths']) == np.ndarray, "Params must include 'stage_lengths' as numpy array'"
                assert type(params['stage_tendons']) == np.ndarray, "Params must include 'stage_tendons' as numpy array"
                assert params['routing'].size()[1] == params['num_tendons'], "Number of stage tendons must equal to number of tendons."
                assert np.size(params['stage_tendons']) == np.size(params['stage_lengths'])*params['num_tendons'], "stage_tendons must have the same size as num_tendons."

            self._valid_build = True
            print("Valid build.")

        except:

            raise ValueError(
                "Building failed! Params must include a valid routing and robot type or include appropriate stage details for multistage robots.")

    # def createEvenlySpaced(self, angle, num_tendons, radius):
    #     """Overwrites the current tendon routing or create tendon routing if it is not in dictionary.
    #         CAUTION: Works for single stage for now only!
    #         To do: do one for multistage.
    #         angle: float (radians)
    #         num_tendons: int 
    #         radius: float (m)

    #    y  ____
    #           |         z points toward into the screen.
    #           |
    #          x

    #     """

    #     assert type(num_tendons) == int, "num_tendons must be an integer!"
        
    #     r = SX.zeros(3, num_tendons)
    #     spacing = 0

    #     for i in range(num_tendons): 

    #         r[0, i] = radius*cos(spacing)
    #         r[1, i] = radius*sin(spacing)
    #         spacing += angle

    #     self.params['num_tendons'] = num_tendons
    #     self.params['routing'] = r

class Robot_Parameter_Builder:

    def __init__(self): 

        self._valid_build = False

    def createCustomParameters(self, params):

        if not self._valid_build: 

            print("Building Continuum Robot using Custom Parameters...")

            self.params = params
            self.physical_type = 'Custom'
            self.params['physical_type'] = self.physical_type

            self._checkDict(self.params)

        else: 

            raise ValueError("Parameters have been built before.")

        return self.params

    def createHollowRod(self, params):


        if not self._valid_build: 

            print("Building Continuum Robot with Hollow Rod...")

            self.params = params
            self._ir = params['inner_radius']
            self._or = params['outer_radius']
            self._E = params['elastic_modulus']
            self._shear_mod = params['shear_modulus']
            try:
                self._rho = params['rho']
                self.mass_distribution = self._rho * self._area
            except:
                self.mass_distribution = params['mass_distribution']
            self._area = np.pi * self._or**2 - np.pi * self._ir**2
            self._I = ((np.pi * self._or**4) / 4) - \
                ((np.pi * self._ir**4) / 4)
            self._J = 2 * self._I
            self.Kse = diag([self._shear_mod * self._area,
                            self._shear_mod * self._area, self._E * self._area])
            self.Kbt = diag([self._E * self._I, self._E *
                            self._I, self._shear_mod * self._J])
            self.physical_type = 'hollow_rod'
            self.params['physical_type'] = self.physical_type
            self.params['Kse'] = self.Kse
            self.params['Kbt'] = self.Kbt
            self.params['mass_distribution'] = self.mass_distribution

            self._checkDict(self.params)

        else: 

            raise ValueError("Parameters have been built before.")


        return self.params

    def createRod(self, params):

        if not self._valid_build: 

            print("Building Continuum Robot with Solid Rod...")

            self.params = params
            self._r = params['radius']
            self._E = params['elastic_modulus']
            self._shear_mod = params['shear_modulus']
            try:
                self._rho = params['rho']
                self.mass_distribution = self._rho * self._area
            except:
                self.mass_distribution = params['mass_distribution']
            self._area = np.pi * self._r**2
            self._I = (np.pi * self._r**4) / 4
            self._J = 2 * self._I
            self.Kse = diag([self._shear_mod * self._area,
                            self._shear_mod * self._area, self._E * self._area])
            self.Kbt = diag([self._E * self._I, self._E *
                            self._I, self._shear_mod * self._J])
            self.mass_distribution = self._rho * self._area
            self.physical_type = 'Solid Rod'
            self.params['physical_type'] = self.physical_type
            self.params['Kse'] = self.Kse
            self.params['Kbt'] = self.Kbt
            self.params['mass_distribution'] = self.mass_distribution

            self._checkDict(self.params)

        else: 

            raise ValueError("Parameters have been built before.")



        return self.params

    def _checkDict(self, params):

        try:

            assert type(params['mass_distribution']) == float, "mass_distribution must be a float."
            assert type(params['tip_weight']) == float, "tip_weight must be a float." 
            assert type(params['Kse']) == DM, "Kse must be of DM type."
            assert type(params['Kbt']) == DM, "Kbt must be of DM type."
            assert type(params['robot_length']) == float, "robot_length must be a float."
            assert type(params['tip_weight']) == float, "tip_weight must be a float."
            self._valid_build = True
            print("Valid build.")

        except:

            raise ValueError(
                "Building failed! Params must include mass distribution, tip weight, Kse, Kbt, robot length.")

if __name__ == "__main__":

    robot_dict = {}
    robot_dict['type'] = 'hollow_rod'
    robot_dict['outer_radius'] = 0.002
    robot_dict['inner_radius'] = 0.0008
    robot_dict['shear_modulus'] = 80e9
    robot_dict['elastic_modulus'] = 3.4698e9
    robot_dict['mass_distribution'] = 0.274
    robot_dict['angle_spacing'] = 2 * pi / 3
    robot_dict['num_tendons'] = 6
    robot_dict['integration_steps'] = 30
    robot_dict['time_step'] = 0.01
    robot_dict['robot_length'] = 0.2
    robot_dict['tendon_offset1'] = 0.01506
    robot_dict['tendon_offset2'] = 0.01506
    robot_dict['tip_weight'] = 0.0

    builder = Robot_Parameter_Builder()
    print(builder.createHollowRod(robot_dict))

    robot_type = {}
    robot_type['routing'] = SX([1, 2, 3])
    robot_type['num_tendons'] = 3
    robot_type['stage_lengths'] = np.array([0.2, 0.2])
    robot_type['stage_tendons'] = np.array([[1, 0, 0], [0, 1, 1]])
    print(robot_type['routing'])

    builder_type = Robot_Type_Builder()
    print(builder_type.createMultiStageStraight(robot_type))
    print(builder_type.params)
