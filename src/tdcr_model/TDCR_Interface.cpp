#include "TDCR_Interface.hpp" 

#define assertm(exp, msg) assert(((void)msg, exp))

TDCR_Interface::TDCR_Interface(MultistageTDCR_Solver* TDCR_, ControllerInterface* MPC_) 

{ 
    MPC = MPC_;
    TDCR = TDCR_;
    desiredTensions.resize(TDCR->getNumTendons(), 1); 
    desiredTensions.setZero();
    deltaTension.resize(TDCR->getNumTendons(), 1); 
    deltaTension.setZero();
    JacobianEta.resize(7, TDCR->getNumTendons());
    I_mxm.resize(TDCR->getNumTendons(), TDCR->getNumTendons()); 
    I_mxm.setZero();
    scaleLoadCell = 1e-3; 

}

Eigen::MatrixXd TDCR_Interface::getTau()

{

    return TDCR->getRobotStates(0).block(18, 0, TDCR->getNumTendons(), 1);

}

Eigen::MatrixXd TDCR_Interface::getDesiredTensions() 

{

    return desiredTensions;

}

Eigen::MatrixXd TDCR_Interface::saveFullStates(std::string fileName) 

{

    TDCR->saveData(fileName, TDCR->integrateFullStates());

    return TDCR->getFullStates();

}

std::vector<Eigen::Matrix<double, 7, 1>> TDCR_Interface::getPoseWorld() 

{

    return TDCR->getRobotPoseWorld();

}

void TDCR_Interface::setDimensions(double numControlStates, std::vector<Eigen::MatrixXi> CSM, Eigen::MatrixXi stagesControlled_) 

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

    stagesControlled = stagesControlled_;
    controlStatesMatrix = CSM;
    unsigned int CSM_size = CSM.size();

    assertm(stagesControlled.rows() == CSM_size, "stagesControlled must equal to controlStatesMatrix size, check dimensions");
    assertm(stagesControlled.cols() == 1, "stagesControlled must equal to controlStatesMatrix size, check dimensions");

    dimensionsSet = true;

}

bool TDCR_Interface::checkBoundaryConditions() 

{

    for (unsigned int i = 0; i < TDCR->getNumStages(); ++i) 
    
    {

        if (TDCR->getBoundaryConditions(i).norm() > 1e-2) 
        
        {

            printf("Boundary Conditions are being violated, please stop the program.");
            return true; 

        }

    }

    return false;
    

}

Eigen::MatrixXd TDCR_Interface::getInitialCondition() 

{

    return TDCR->getInitialConditions();

}

Eigen::MatrixXd TDCR_Interface::getHighLevelControl(Eigen::MatrixXd poseDesired, Eigen::MatrixXd currentPose) 

{

    /* This method is a generalised version of the tip controller. 
    */ 

    assertm(currentPose.cols() == 1, "currentPose must be of TDCR->getNumTendons() + poseDesired.rows() by 1.");
    assertm(currentPose.rows() == TDCR->getNumTendons() + poseDesired.rows(), "currentPose must be of TDCR->getNumTendons() + poseDesired.rows() by 1.");
    assertm(dimensionsSet == true, "Please set the dimensions required first.");
    assertm(poseDesired.rows() == customPose.rows(), "poseDesired must be of numControlStates by 1. ");
    assertm(poseDesired.cols() == 1, "poseDesired must be of numControlStates by 1. "); 

    unsigned int numStatesAtStage_i; 
    unsigned int numStatesCumulative = 0; 
    unsigned int index; 
    Eigen::Matrix<double, 7, 1> poseAtStage_i;
    Eigen::MatrixXd poseDesiredMPC(TDCR->getNumTendons() + poseDesired.rows(), 1);
    Eigen::MatrixXd TendonsDesiredMPC(TDCR->getNumTendons(), 1);
    unsigned int stage_num;
    TendonsDesiredMPC.setZero();

    for (unsigned int i = 0; i < stagesControlled.rows(); ++i) 
    
    {

        stage_num = stagesControlled(i, 0);
        numStatesAtStage_i = controlStatesMatrix[i].rows();
        poseAtStage_i = TDCR->getRobotPoseWorld(stage_num);

        for (unsigned int k = 0; k < numStatesAtStage_i; ++k) 
        
        {

            index = controlStatesMatrix[i](k, 0);
            customJacobianEta.block(numStatesCumulative + k, 0, 1, TDCR->getNumTendons()) = TDCR->getJacobiansEta()[stage_num].row(index);
            customPose(numStatesCumulative + k, 0) = poseAtStage_i(index, 0);

        }
 
        numStatesCumulative += numStatesAtStage_i;

    }

    poseDesiredMPC << poseDesired, TendonsDesiredMPC;
    customPoseError = poseDesired - customPose; 

    deltaTension = MPC->solveOptimalControl(currentPose, customJacobianEta, poseDesiredMPC);

    desiredTensions += deltaTension*TDCR->getSamplingTime();

    return deltaTension;


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
    unsigned int stage_num;
    Eigen::Matrix<double, 7, 1> poseAtStage_i;
    Eigen::MatrixXd poseDesiredMPC(TDCR->getNumTendons() + poseDesired.rows(), 1);
    Eigen::MatrixXd TendonsDesiredMPC(TDCR->getNumTendons(), 1);
    Eigen::MatrixXd customPoseMPC(TDCR->getNumTendons() + poseDesired.rows(), 1);
    TendonsDesiredMPC.setZero();

    for (unsigned int i = 0; i < stagesControlled.rows(); ++i) 
    
    {

        // i is the stage number. 

        stage_num = stagesControlled(i, 0);
        numStatesAtStage_i = controlStatesMatrix[i].rows();
        poseAtStage_i = TDCR->getRobotPoseWorld(stage_num);

        for (unsigned int k = 0; k < numStatesAtStage_i; ++k) 
        
        {

            index = controlStatesMatrix[i](k, 0);
            customJacobianEta.block(numStatesCumulative + k, 0, 1, TDCR->getNumTendons()) = TDCR->getJacobiansEta()[stage_num].row(index);
            customPose(numStatesCumulative + k, 0) = poseAtStage_i(index, 0);

        }
 
        numStatesCumulative += numStatesAtStage_i;

    }

    poseDesiredMPC << poseDesired, TendonsDesiredMPC;
    customPoseMPC << customPose, TDCR->getTau();
    customPoseError = poseDesired - customPose; 

    deltaTension = MPC->solveOptimalControl(customPoseMPC, customJacobianEta, poseDesiredMPC);

    desiredTensions += deltaTension*TDCR->getSamplingTime();

    return deltaTension;

}

void TDCR_Interface::setScaleLoadCell(double value)

{

    scaleLoadCell = value; 

}

Eigen::MatrixXd TDCR_Interface::getCustomPoseError()

{

    return customPoseError;

}

Eigen::MatrixXd TDCR_Interface::getCustomPose()

{

    return customPose;

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