#include "TDCR_Interface.hpp" 

#define assertm(exp, msg) assert(((void)msg, exp))

TDCR_Interface::TDCR_Interface(MultistageTDCR_Solver& TDCR_) 

{ 

    TDCR = &TDCR_;
    desiredTensions.resize(TDCR->getNumTendons(), 1); 
    desiredTensions.setZero();
    deltaTension.resize(TDCR->getNumTendons(), 1); 
    deltaTension.setZero();
    weightPriorityTip.resize(7, 7);
    weightPriorityTip.setIdentity();
    JacobianEta.resize(7, TDCR->getNumTendons());
    I_mxm.resize(TDCR->getNumTendons(), TDCR->getNumTendons()); 
    I_mxm.setZero();
    lambda = 0.1; 
    Kp = 1.; 
    scaleLoadCell = 1e-3; 

}

Eigen::MatrixXd TDCR_Interface::getDesiredTensions() 

{

    return desiredTensions;

}

void TDCR_Interface::setDimensions(double numControlStates, std::vector<Eigen::MatrixXi> CSM) 

{

    // Check if number of elements in CSM = numControlStates.

    int num_elements = 0; 

    for (unsigned int i = 0; i < TDCR->getNumStages(); ++i) 
    
    {

        num_elements += CSM[i].rows(); 

    }

    assertm(num_elements == numControlStates, "Number of elements in CSM must be the same as numControlStates.");

    customJacobianEta.resize(numControlStates, TDCR->getNumTendons());
    customPose.resize(numControlStates, 1);
    customPoseError.resize(numControlStates, 1);
    weightPriorityCustom.resize(numControlStates, numControlStates);
    weightPriorityCustom.setIdentity();

    controlStatesMatrix = CSM;
    dimensionsSet = true;

}

bool TDCR_Interface::checkBoundaryConditions() 

{

    for (unsigned int i = 0; i < TDCR->getNumStages(); ++i) 
    
    {

        if (TDCR->getBoundaryConditions(i).norm() < 1e-2) 
        
        {

            return true; 

        }

    }

    return false;
    

}

Eigen::MatrixXd TDCR_Interface::getHighLevelControl(Eigen::MatrixXd poseDesired) 

{

    /* This method is a generalised version of the tip controller. 
    */ 

    assertm(dimensionsSet == true, "Please set the dimensions required first.");
    assertm(poseDesired.rows() == customPose.rows(), "poseDesired must be of numControlStates by 1. ");
    assertm(poseDesired.cols() == 1, "poseDesired must be of numControlStates by 1. "); 

    unsigned int numStatesAtStage_i; 
    unsigned int numStatesCumulative = 0; 
    unsigned int index; 
    Eigen::Matrix<double, 7, 1> poseAtStage_i;

    for (unsigned int i = 0; i < TDCR->getNumStages(); ++i) 
    
    {

        numStatesAtStage_i = controlStatesMatrix[i].rows();
        poseAtStage_i = MathUtils::robotStates2Pose(TDCR->getRobotStates(i).block<12, 1>(0, 0));

        for (unsigned int k = 0; k < numStatesAtStage_i; ++k) 
        
        {

            index = controlStatesMatrix[i](k, 0);
            customJacobianEta.block(numStatesCumulative + k, 0, 1, TDCR->getNumTendons()) = TDCR->getJacobiansEta()[i].row(index);
            customPose(numStatesCumulative + k, 0) = poseAtStage_i(index, 0);

        }
 
        numStatesCumulative += numStatesAtStage_i;

    }

    customPoseError = poseDesired - customPose; 

    deltaTension = Kp*(customJacobianEta.transpose() * weightPriorityCustom * customJacobianEta + pow(lambda, 2)*I_mxm).completeOrthogonalDecomposition().pseudoInverse() * customJacobianEta.transpose() * weightPriorityCustom * customPoseError;

    std::cout << "deltaTension: " << (customJacobianEta.transpose() * weightPriorityCustom * customJacobianEta + pow(lambda, 2)*I_mxm).inverse()  << "\n\n";
    std::cout << weightPriorityCustom << "\n\n";

    desiredTensions += deltaTension*TDCR->getSamplingTime();

    return desiredTensions;

}

Eigen::MatrixXd TDCR_Interface::getHighLevelControl(Eigen::Matrix<double, 7, 1> poseDesired) 

{

    JacobianEta = TDCR->getJacobiansEta()[TDCR->getNumStages()];

    poseTip = MathUtils::robotStates2Pose(TDCR->getRobotStates(TDCR->getNumStages()));
    poseTipError = poseDesired - poseTip;

    deltaTension = Kp*(JacobianEta.transpose() * weightPriorityTip * JacobianEta + pow(lambda, 2)*I_mxm).inverse() * JacobianEta.transpose() * weightPriorityTip * poseTipError ;

    desiredTensions += deltaTension*TDCR->getSamplingTime();

    return desiredTensions;

}

void TDCR_Interface::setScaleLoadCell(double value)

{

    scaleLoadCell = value; 

}

void TDCR_Interface::setWeightsAllStages(Eigen::MatrixXd W_custom)

{

    assertm(W_custom.rows() == weightPriorityCustom.rows(), "W_custom must be of size 7 by num_tendons");
    assertm(W_custom.cols() == weightPriorityCustom.cols(), "W_custom must be of size 7 by num_tendons");
    weightPriorityCustom = W_custom;

}

void TDCR_Interface::setWeightsTip(Eigen::MatrixXd W) 

{

    assertm(W.rows() == weightPriorityTip.rows(), "W must be of size 7 by num_tendons");
    assertm(W.cols() == weightPriorityTip.cols(), "W must be of size 7 by num_tendons");
    weightPriorityTip = W; 

}

Eigen::MatrixXd TDCR_Interface::getCustomPoseError()

{

    return customPoseError;

}

Eigen::Matrix<double, 7, 1> TDCR_Interface::getPoseError()

{

    return poseTipError;

}


void TDCR_Interface::simulateStep(Eigen::MatrixXd tau) 

{

    TDCR->simulateStep(tau);

}


void TDCR_Interface::trackMeasuredTension(Eigen::MatrixXd tauTrack) 

{

    assertm(tauTrack.rows() == TDCR->getNumTendons(), "tauTrack must be TDCR->getNumTendons() by 1.");
    assertm(tauTrack.cols() == 1, "tauTrack must be TDCR->getNumTendons() by 1.");


    Eigen::MatrixXd tauError(TDCR->getNumTendons(), 1);
    Eigen::MatrixXd inputStep(TDCR->getNumTendons(), 1); 

    tauError = tauTrack - TDCR->getTau();
    inputStep = tauError * 9.81 * scaleLoadCell / TDCR->getSamplingTime();


    TDCR->simulateStep(inputStep);


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

    tau.setZero();
    TDCR->simulateStep(tau);

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

std::vector<Eigen::MatrixXd> TDCR_Interface::getJacobiansEta() 

{

    return TDCR->getJacobiansEta();

}


TDCR_Interface::~TDCR_Interface() {}; 