#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <IntegratorInterface.hpp>
#include <vector>
#include <cassert>

#define assertm(exp, msg) assert(((void)msg, exp))
#define EPS 1e-9;

class MultistageTDCR_Solver {

    public: 

        MultistageTDCR_Solver(int numTendons, int numStages, std::vector<IntegrationInterface> integrators, std::vector<IntegrationInterface> integratorsStep, Eigen::MatrixXd stage_tendons, Eigen::MatrixXd routing_);
        virtual ~MultistageTDCR_Solver();
        std::vector<Eigen::MatrixXd> getRobotStates(bool print_level, bool csv_level);
        void solveForwardKinematics();
        void simulateStep();
        void setTau(Eigen::MatrixXd tau);

    private: 

        void setNumTendons(int num, Eigen::MatrixXd routing_);
        void setNumStages(int num);
        void initialiseJacobianMatrices(Eigen::MatrixXd stage_tendons_);
        void setInitialConditions(unsigned int stage_num, Eigen::Matrix<double, 6, 1> ic_force_moment);
        void forwardFiniteDifferences(unsigned int stage_num);  
        void convertStageTendonsIndex();

        std::vector<Eigen::MatrixXd> robotStates;
        std::vector<Eigen::MatrixXd> initialConditions;
        std::vector<IntegrationInterface> integrators;
        std::vector<IntegrationInterface> integratorsStep;
        std::vector<Eigen::VectorXi> stageTendonsIndex;
        std::vector<Eigen::MatrixXd> J_q, E_q, B_q, B_yu, E_yu, B;

        Eigen::MatrixXd integrateStates(unsigned int stage_num); 
        Eigen::MatrixXd integrateStep(unsigned int stage_num); 


        Eigen::MatrixXd solveJacobians();
        Eigen::Matrix<double, 6, 1> getBoundaryConditions(unsigned int stage_num);
        Eigen::Matrix<double, 6, 1> getPointForceMoment(unsigned int stage_num);
        Eigen::MatrixXd integrateWithIncrement(unsigned int index, unsigned int stage_num);

        Eigen::MatrixXd stageTendons;
        Eigen::Matrix<double, 3, 1> v; 
        Eigen::MatrixXd R; 
        Eigen::MatrixXd x;
        Eigen::MatrixXd routing; 
        Eigen::VectorXd stageLengths;
        
        const static int num_p = 3; 
        const static int num_R = 9; 
        const static int num_m = 3; 
        const static int num_n = 3; 
        unsigned int num_tendons; 
        const static int num_alpha = 1; 
        unsigned int num_stages;
        unsigned int num_total;

};