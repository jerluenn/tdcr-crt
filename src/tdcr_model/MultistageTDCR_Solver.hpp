#ifndef MULTISTAGETDCR_H
#define MULTISTAGETDCR_H

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <IntegratorInterface.hpp>
#include <MathUtils.hpp>
#include <vector>
#include <Integrator.hpp>
#include <cassert>
#include <chrono>



class MultistageTDCR_Solver {

    public: 

        MultistageTDCR_Solver(int numTendons, int numStages, std::vector<IntegrationInterface> integrators, std::vector<IntegrationInterface> integratorsStep, Eigen::MatrixXd stage_tendons, Eigen::MatrixXd routing_);
        MultistageTDCR_Solver();
        virtual ~MultistageTDCR_Solver();
        std::vector<Eigen::MatrixXd> getRobotStates(bool print_level);
        Eigen::MatrixXd getRobotStates(unsigned int stage_num);
        void simulateStep(Eigen::MatrixXd tau);
        void setTau(Eigen::MatrixXd tau);
        void setInitialConditions(unsigned int stage_num, Eigen::Matrix<double, 6, 1> ic_force_moment);
        unsigned int getNumStages(); 
        unsigned int getNumTendons();
        Eigen::MatrixXd getTau(); 
        void integrateStatesandUpdate(unsigned int stage_num);
        std::vector<Eigen::MatrixXd> getJacobians();
        std::vector<Eigen::MatrixXd> getJacobiansEta();
        double getSamplingTime();
        Eigen::MatrixXd integrateStates(unsigned int stage_num);
        Eigen::MatrixXd integrateStep(unsigned int stage_num); 
        Eigen::Matrix<double, 6, 1> getBoundaryConditions(unsigned int stage_num, Eigen::MatrixXd robotStates_, Eigen::MatrixXd initialConditions_);
        Eigen::Matrix<double, 6, 1> getBoundaryConditions(unsigned int stage_num, Eigen::MatrixXd robotStates_);
        Eigen::Matrix<double, 6, 1> getBoundaryConditions(unsigned int stage_num);

        MathUtils::Timer timer;
        

    private: 

        void setNumTendons(int num, Eigen::MatrixXd routing_);
        void setNumStages(int num);
        void initialiseJacobianMatrices(Eigen::MatrixXd stage_tendons_);
        void convertStageTendonsIndex();
        void testFunction();

        std::vector<Eigen::MatrixXd> robotStates;
        std::vector<Eigen::MatrixXd> initialConditions;
        std::vector<IntegrationInterface> integrators;
        std::vector<IntegrationInterface> integratorsStep;
        std::vector<Eigen::VectorXi> stageTendonsIndex_Termination;
        std::vector<Eigen::VectorXi> stageTendonsIndex_Jacobians;
        std::vector<Eigen::MatrixXd> J_q, E_q, B_q, B_yu, E_yu, J_b, E_yu_Nplus1, B_yu_Nplus1;
        std::vector<Eigen::MatrixXd> J_world, J_world_eta;

        void solveJacobians();

        Eigen::Matrix<double, 6, 1> getPointForceMoment(unsigned int stage_num);
        Eigen::Matrix<double, 6, 1> getPointForceMoment(unsigned int stage_num, Eigen::MatrixXd robotStates_);
        Eigen::Affine3d getAffine3d_Distal(unsigned int stage_num);
        Eigen::Affine3d getAffine3d_Proximal(unsigned int stage_num);
        Eigen::MatrixXd integrateWithIncrement(unsigned int index, unsigned int stage_num);
        Eigen::MatrixXd addIncrement(unsigned int index, unsigned int stage_num);
        Eigen::MatrixXd convertJacobianToJacobianEta(unsigned int index, Eigen::MatrixXd J_world_tmp, Eigen::Matrix3d R_tmp);

        Eigen::Matrix3d I_3x3;
        Eigen::MatrixXd stageTendons;
        Eigen::Matrix<double, 6, 6> AdjointMatrix;
        Eigen::Matrix<double, 3, 1> v; 
        Eigen::Affine3d T_Affine3d; 
        Eigen::MatrixXd R, eta_dot; 
        Eigen::MatrixXd x; // dummy variable to hold all the states at initial condition.
        Eigen::MatrixXd routing; 
        Eigen::VectorXd stageLengths;
        const static int num_p = 3; 
        const static int num_R = 9; 
        const static int num_m = 3; 
        const static int num_n = 3; 
        unsigned int num_tendons; 
        constexpr static double EPS = 1.0e-9;
        double dt = 0.005;
        const static int num_alpha = 1; 
        unsigned int num_stages;
        unsigned int num_total;

};

#endif