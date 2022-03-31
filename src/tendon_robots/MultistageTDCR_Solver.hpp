#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <Eigen/Dense>
#include <IntegratorInterface.hpp>
#include <vector>
#include <cassert>

#define assertm(exp, msg) assert(((void)msg, exp))

class MultistageTDCR_Solver {


    public: 

        MultistageTDCR_Solver(int numTendons, int numStages, std::vector<IntegrationInterface> integrators, std::vector<IntegrationInterface> integratorsStep);
        virtual ~MultistageTDCR_Solver();
        std::vector<Eigen::MatrixXd> getRobotStates(bool print_level, bool csv_level);
        void solveForwardKinematics();
        void simulateStep();

    private: 

        void setNumTendons(int num);
        void setNumStages(int num);
        std::vector<Eigen::MatrixXd> robotStates;
        std::vector<Eigen::MatrixXd> initialConditions;
        std::vector<IntegrationInterface> integrators;
        std::vector<IntegrationInterface> integratorsStep;
        Eigen::MatrixXd integrateStates(unsigned int stage_num); 
        Eigen::MatrixXd integrateStep(unsigned int stage_num); 
        void setInitialConditions(unsigned int stage_num, Eigen::Matrix<double, 6, 1> ic_force_moment);
        Eigen::MatrixXd x;
        const static int num_p = 3; 
        const static int num_R = 9; 
        const static int num_m = 3; 
        const static int num_n = 3; 
        const static int num_alpha = 1; 
        unsigned int num_stages;
        unsigned int num_tendons; 
        unsigned int num_total;

};