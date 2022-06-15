#include "TDCR_Interface.hpp" 

TDCR_Interface::TDCR_Interface(MultistageTDCR_Solver& TDCR_) 

{ 

    TDCR = &TDCR_;
    desiredTension.resize(TDCR->getNumTendons()); 
    desiredTension.setZero();
    deltaTension.resize(TDCR->getNumTendons()); 
    deltaTension.setZero();

}

Eigen::MatrixXd TDCR_Interface::getHighLevelControl(Eigen::Matrix<double, 7, 1> desiredPose) 

{

    Eigen::Matrix<double, 7, 1> pose;
    Eigen::Matrix<double, 7, 1> errorPose;

    pose = MathUtils::robotStates2Pose(TDCR->getRobotStates(TDCR->getNumStages()));
    errorPose = desiredPose - pose;

    // delta_q = 5*JacobianEta.transpose() * (JacobianEta * JacobianEta.transpose() 
    // + pow(lambda, 2)*I_7x7).inverse() * poseError;

    // Need to first get JacobianEta from MultistageTDCR side..

    deltaTension += delta_q;

    // Update desired tensions to send to motor PID.  

}

Eigen::MatrixXd TDCR_Interface::getHighLevelControl(std::vector<Eigen::Matrix<double, 7, 1>> desiredPose) 

{

    

}

void TDCR_Interface::trackMeasuredTension(Eigen::MatrixXd tau) 

{

    // Assert tau size. 

    // Get tauError. 

    // Convert to Newtons.

    // Input to simulation = tauError/timeStep (Write getTimeStep() fn.)

    // SimulateStep()
    

}

void TDCR_Interface::solveForwardKinematics(Eigen::MatrixXd tau, bool print_level) 

{

    assertm(tau.rows() == TDCR->getNumTendons(), "tau must have num_tendons elements"); 
    assertm(tau.cols() == 1, "tau must be a vector!");
    TDCR->setTau(tau);

    LevenbergMarquardtFunctor LMFunctor(*TDCR);
    Eigen::VectorXd x(TDCR->getNumStages()*6);
    x.setZero();

	Eigen::LevenbergMarquardt<LevenbergMarquardtFunctor, double> lm(LMFunctor);
    lm.minimizeInit(x);
	int status = lm.minimize(x);

    for (unsigned int i = 0; i < TDCR->getNumStages(); ++i) 
    
    {

        TDCR->setInitialConditions(i, x.segment<6>(i*6));
        TDCR->integrateStatesandUpdate(i);

    }

    if (print_level == true) 

    {

        std::cout << "status: " << status << "\n\n";
        std::cout << "fnorm: " << lm.fnorm << "\n\n";
        std::cout << "Sol: " << x << std::endl;

    }


}

std::vector<Eigen::MatrixXd> TDCR_Interface::getJacobians() 

{

    return TDCR->getJacobians();

}

TDCR_Interface::~TDCR_Interface() {}; 