#include "MultistageTDCR_Solver.hpp"

Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

MultistageTDCR_Solver::MultistageTDCR_Solver(int numTendons, int numStages, std::vector<IntegrationInterface> integrators_, std::vector<IntegrationInterface> integratorsStep_)
{               

    integrators = integrators_;
    integratorsStep = integratorsStep_;
    setNumTendons(numTendons);
    setNumStages(numStages);
    initialiseJacobianMatrices();


    Eigen::Matrix<double, 6, 1> ic; 
    ic << 5.91802264e-01, -2.76353034e-29, -5.00264992e+00,  8.99753480e-29,
  1.44320521e-01, -4.97943595e-29;
    setInitialConditions(0, ic);

    robotStates[0] = integrateStates(0);


}

void MultistageTDCR_Solver::convertStageTendonsIndex() 

{}

void MultistageTDCR_Solver::initialiseJacobianMatrices() {

    Eigen::MatrixXd tmp; 
    unsigned int tmp_int = num_tendons;

    for(int i = 0; i < num_stages; ++i) {

        tmp.resize(6,tmp_int);
        J_q.push_back(tmp);
        tmp.resize(6,tmp_int);
        E_q.push_back(tmp);
        tmp.resize(6,tmp_int);
        B_q.push_back(tmp);
        tmp.resize(6,6);
        B_yu.push_back(tmp);
        tmp.resize(6,6);
        E_yu.push_back(tmp);

        tmp_int -= stage_tendons.row(i).sum();

    }



}

void MultistageTDCR_Solver::forwardFiniteDifferences(unsigned int stage_num) {



}

Eigen::MatrixXd MultistageTDCR_Solver::integrateWithIncrement(unsigned int index, unsigned int stage_num) {

    x = initialConditions[stage_num];
    x(index) += EPS;
    return integrators[stage_num].integrate(x);

}

Eigen::MatrixXd MultistageTDCR_Solver::getBoundaryConditions(unsigned int stage_num, Eigen::MatrixXd integrated_states) 

{



} 

Eigen::MatrixXd MultistageTDCR_Solver::getPointForceMoment(unsigned int stage_num, Eigen::MatrixXd integrated_states) 

{



} 


void MultistageTDCR_Solver::setTau(Eigen::MatrixXd tau) 

{

    assertm(num_tendons == tau.rows(), "tau must be the same size as the number of stages!");

    for (unsigned int i = 0; i < num_stages; ++i) {

        initialConditions[i].block(num_p + num_R + num_m + num_n, 0, num_tendons, 1) = tau;

    }

}


void MultistageTDCR_Solver::setNumTendons(int num) 

{

    num_tendons = num;
    num_total = num_p + num_R + num_m + num_n + num_alpha + num_tendons;

    x.resize(num_total, 1);

}

void MultistageTDCR_Solver::setNumStages(int num) 

{

    num_stages = num;
    Eigen::MatrixXd states;
    states.resize(num_total, 1);
    states.setZero();
    states.data()[18] = 5.0;
    states.block<9, 1>(3, 0) << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;

    for (unsigned int i = 0; i < num_stages; ++i ) {

        initialConditions.push_back(states);
        robotStates.push_back(states);

    } 


    assertm(num_stages == integrators.size(), "Checking integrator list size."); 
    assertm(num_stages == integratorsStep.size(), "Checking integrator list size.");

}


std::vector<Eigen::MatrixXd> MultistageTDCR_Solver::getRobotStates(bool print_level, bool csv_level)
{


    if (print_level == true) {

        for (unsigned int i = 0; i < num_stages; ++i) {

            std::cout << "Initial conditions at stage " << i << ": \n";
            std::cout << initialConditions[i].transpose().format(OctaveFmt) << "\n";
            std::cout << "Robot states at stage " << i << ": \n";
            std::cout << robotStates[i].transpose().format(OctaveFmt) << "\n";

        }

    } 

    return robotStates;

}

Eigen::MatrixXd MultistageTDCR_Solver::integrateStates(unsigned int stage_num)
{

    return integrators[stage_num].integrate(initialConditions[stage_num]);

} 

Eigen::MatrixXd MultistageTDCR_Solver::integrateStep(unsigned int stage_num)

{

    return integrators[stage_num].integrate(initialConditions[stage_num]);

}

void MultistageTDCR_Solver::solveForwardKinematics()
{


}

void MultistageTDCR_Solver::setInitialConditions(unsigned int stage_num, Eigen::Matrix<double, 6, 1> ic_force_moment) {

    initialConditions[stage_num].block<6,1>(num_p + num_R, 0) = ic_force_moment;

} 

MultistageTDCR_Solver::~MultistageTDCR_Solver() {};

