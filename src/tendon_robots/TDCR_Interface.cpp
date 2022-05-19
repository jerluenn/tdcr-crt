#include "TDCR_Interface.hpp" 

TDCR_Interface::TDCR_Interface(MultistageTDCR_Solver& TDCR_) 

{ 

    
    TDCR = &TDCR_;

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

TDCR_Interface::~TDCR_Interface() {}; 