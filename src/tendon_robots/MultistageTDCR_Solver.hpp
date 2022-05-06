#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <IntegratorInterface.hpp>
#include <MathUtils.hpp>
#include <vector>
#include <Integrator.hpp>
#include <cassert>
#include <chrono>

#define assertm(exp, msg) assert(((void)msg, exp))

class MultistageTDCR_Solver {

    public: 

        MultistageTDCR_Solver(int numTendons, int numStages, std::vector<IntegrationInterface> integrators, std::vector<IntegrationInterface> integratorsStep, Eigen::MatrixXd stage_tendons, Eigen::MatrixXd routing_);
        virtual ~MultistageTDCR_Solver();
        std::vector<Eigen::MatrixXd> getRobotStates(bool print_level, bool csv_level);
        void solveForwardKinematics();
        void simulateStep(Eigen::MatrixXd tau);
        void setTau(Eigen::MatrixXd tau);
        

    private: 

        MathUtils::Timer timer;

        void setNumTendons(int num, Eigen::MatrixXd routing_);
        void setNumStages(int num);
        void initialiseJacobianMatrices(Eigen::MatrixXd stage_tendons_);
        void setInitialConditions(unsigned int stage_num, Eigen::Matrix<double, 6, 1> ic_force_moment);
        void convertStageTendonsIndex();

        std::vector<Eigen::MatrixXd> robotStates;
        std::vector<Eigen::MatrixXd> initialConditions;
        std::vector<IntegrationInterface> integrators;
        std::vector<IntegrationInterface> integratorsStep;
        std::vector<Eigen::VectorXi> stageTendonsIndex_Termination;
        std::vector<Eigen::VectorXi> stageTendonsIndex_Jacobians;
        std::vector<Eigen::MatrixXd> J_q, E_q, B_q, B_yu, E_yu, J_b, E_yu_Nplus1, B_yu_Nplus1;
        std::vector<Eigen::MatrixXd> J_world;

        Eigen::MatrixXd integrateStates(unsigned int stage_num); 
        Eigen::MatrixXd integrateStep(unsigned int stage_num); 

        void solveJacobians();
        Eigen::Matrix<double, 6, 1> getBoundaryConditions(unsigned int stage_num, Eigen::MatrixXd robotStates_, Eigen::MatrixXd initialConditions_);
        Eigen::Matrix<double, 6, 1> getBoundaryConditions(unsigned int stage_num, Eigen::MatrixXd robotStates_);
        Eigen::Matrix<double, 6, 1> getBoundaryConditions(unsigned int stage_num);
        Eigen::Matrix<double, 6, 1> getPointForceMoment(unsigned int stage_num);
        Eigen::Matrix<double, 6, 1> getPointForceMoment(unsigned int stage_num, Eigen::MatrixXd robotStates_);
        Eigen::MatrixXd integrateWithIncrement(unsigned int index, unsigned int stage_num);
        Eigen::MatrixXd addIncrement(unsigned int index, unsigned int stage_num);

        Eigen::MatrixXd stageTendons;
        Eigen::Matrix<double, 3, 1> v; 
        Eigen::MatrixXd R; 
        Eigen::MatrixXd x; // dummy variable to hold all the states at initial condition.
        Eigen::MatrixXd routing; 
        Eigen::VectorXd stageLengths;
        const static int num_p = 3; 
        const static int num_R = 9; 
        const static int num_m = 3; 
        const static int num_n = 3; 
        unsigned int num_tendons; 
        constexpr static double EPS = 1.0e-10;
        constexpr static double dt = 0.005;
        const static int num_alpha = 1; 
        unsigned int num_stages;
        unsigned int num_total;

};