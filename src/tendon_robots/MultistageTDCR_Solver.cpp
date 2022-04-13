#include "MultistageTDCR_Solver.hpp"

Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

MultistageTDCR_Solver::MultistageTDCR_Solver(int numTendons, int numStages, std::vector<IntegrationInterface> integrators_, std::vector<IntegrationInterface> integratorsStep_, Eigen::MatrixXd stage_tendons, Eigen::MatrixXd routing_)
{               

    v << 0.0, 0.0, 1.0;
    integrators = integrators_;
    integratorsStep = integratorsStep_;
    setNumTendons(numTendons, routing_);
    setNumStages(numStages);
    initialiseJacobianMatrices(stage_tendons);


    Eigen::Matrix<double, 6, 1> ic; 
    ic << -1.2628939, 1.491191, -4.60228383, 0.0646629, 0.03807737, -0.00540438;
    setInitialConditions(0, ic);
    Eigen::Matrix<double, 6, 1> ic1;
    ic1 << -2.93146117,  1.41209708, -3.83800553,  0.05841705,  0.09182832, -0.01998003;
    setInitialConditions(1, ic1);
    Eigen::MatrixXd tau;
    tau.resize(6, 1); 
    tau << 0, 0, 0, 0, 5, 0;
    setTau(tau);

    robotStates[0] = integrateStates(0);
    robotStates[1] = integrateStates(1);

    unsigned int n = 0;

    std::cout << getBoundaryConditions(n) << "\n\n\n";
    std::cout << getPointForceMoment(n) << "\n\n\n";


}

void MultistageTDCR_Solver::convertStageTendonsIndex() 

{

    std::vector<int> indices; 
    Eigen::VectorXi indices_eig; 

    for (auto it : stageTendons.rowwise()) 

        {

        for (int i = 0; i < it.size();  ++i) {

            if (it(i) == 1) {

                indices.push_back(i);

            } 

        
        }

        indices_eig = Eigen::VectorXi::Map(indices.data(), indices.size());
        stageTendonsIndex.push_back(indices_eig);
        indices.clear();
        
        }


}

void MultistageTDCR_Solver::initialiseJacobianMatrices(Eigen::MatrixXd stage_tendons) {

    assertm(stage_tendons.rows() == num_stages, "Number of rows in stage tendons must equal number of stages!");
    assertm(stage_tendons.cols() == num_tendons, "Number of cols in stage tendons must equal number of tendons!");
    stageTendons = stage_tendons;    
    Eigen::MatrixXd tmp; 
    unsigned int tmp_int = num_tendons;

    convertStageTendonsIndex();

    for(unsigned int i = 0; i < num_stages; ++i) {

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

Eigen::Matrix<double, 6, 1> MultistageTDCR_Solver::getBoundaryConditions(unsigned int stage_num) 

{

    Eigen::Matrix<double, 6, 1> boundaryConditions;
    Eigen::Matrix<double, 6, 1> pointForceMoment;
    Eigen::Matrix<double, 6, 1> internalForcesandMoments_N_0; // internal wrench at stage N, s = 0.
    Eigen::Matrix<double, 6, 1> internalForcesandMoments_Nminus1_l; // internal wrench at stage N-1, s = l. 
    Eigen::Matrix<double, 6, 6> R_diag; 
    Eigen::Matrix<double, 3, 3> zeros_3x3; 
    zeros_3x3.setZero();
    R.resize(9,1);
    pointForceMoment.setZero();

    R = robotStates[stage_num].block<9,1>(3,0); 
    R.resize(3,3); 
    internalForcesandMoments_N_0 = initialConditions[stage_num].block<6, 1>(num_p + num_R, 0);
    internalForcesandMoments_Nminus1_l = robotStates[stage_num].block<6, 1>(num_p + num_R, 0);
    pointForceMoment = getPointForceMoment(stage_num);
    R_diag << R, zeros_3x3, zeros_3x3, R; 

    boundaryConditions = R_diag * internalForcesandMoments_N_0 - internalForcesandMoments_Nminus1_l - pointForceMoment;

    return boundaryConditions;

} 

Eigen::Matrix<double, 6, 1> MultistageTDCR_Solver::getPointForceMoment(unsigned int stage_num) 

{

    /* Get Point Force and Point Moment in the local frame. (Frame N) */

    Eigen::Matrix<double, 6, 1> pointForceMoment;
    R.resize(9,1);
    pointForceMoment.setZero();
    Eigen::Vector3d PointForce;
    Eigen::Vector3d tmpVector; 


    for (auto k : stageTendonsIndex[stage_num]) 
    {

        R = robotStates[stage_num].block<9,1>(3,0); 
        R.resize(3,3);
        PointForce = R*v*initialConditions[stage_num].block(num_p + num_R + num_m + num_n + k, 0, 1, 1);
        pointForceMoment.block<3, 1>(0, 0) -= PointForce;// Point Force.
        tmpVector = (R*routing.col(k));
        pointForceMoment.block<3, 1>(3, 0) -= tmpVector.cross(PointForce); // Point Moment
        
    }
              
    return pointForceMoment;

} 

void MultistageTDCR_Solver::setTau(Eigen::MatrixXd tau) 

{

    assertm(num_tendons == tau.rows(), "tau must be the same size as the number of stages!");

    for (unsigned int i = 0; i < num_stages; ++i) {

        initialConditions[i].block(num_p + num_R + num_m + num_n, 0, num_tendons, 1) = tau;

    }

}


void MultistageTDCR_Solver::setNumTendons(int num, Eigen::MatrixXd routing_) 

{

    num_tendons = num;
    num_total = num_p + num_R + num_m + num_n + num_alpha + num_tendons;

    x.resize(num_total, 1);
    routing.resize(3, num_total);
    routing = routing_;

}

void MultistageTDCR_Solver::setNumStages(int num) 

{

    num_stages = num;
    Eigen::MatrixXd states;
    states.resize(num_total, 1);
    states.setZero();
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

