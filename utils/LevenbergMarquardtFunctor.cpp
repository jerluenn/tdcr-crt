#include "LevenbergMarquardtFunctor.hpp"
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>

LevenbergMarquardtFunctor::LevenbergMarquardtFunctor(MultistageTDCR_Solver& TDCR_, unsigned int stage_num_)

{

    TDCR = &TDCR_;
    m = 6; 
    n = 6;
    stage_num = stage_num_;
    // Insert stage_num.
    // Insert vector x at the distal end of N. 
    // Update boundary conditions.


}

LevenbergMarquardtFunctor::LevenbergMarquardtFunctor() 

{}

int LevenbergMarquardtFunctor::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const 

{

    TDCR->setInitialConditions(stage_num, x.segment<6>(0));
    TDCR->integrateStatesandUpdate(stage_num); 
    fvec.segment<6>(0) = TDCR->getBoundaryConditions(stage_num);

    return 0;

}

int LevenbergMarquardtFunctor::df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const 

{

    Eigen::MatrixXd f, f_plus; 
    

    f.resize(6, 1);
    f_plus.resize(6, 1);

    TDCR->setInitialConditions(stage_num, x.segment<6>(0));
    TDCR->integrateStatesandUpdate(stage_num); 
    f.block<6, 1>(0, 0) = TDCR->getBoundaryConditions(stage_num); 

    for (unsigned int i = 0; i < 6; ++i) 

    {

        Eigen::VectorXd x_plus(x);
        x_plus(i) += EPS;

        TDCR->setInitialConditions(stage_num, x_plus.segment<6>(0));
        TDCR->integrateStatesandUpdate(stage_num);           
        f_plus.block<6, 1>(0, 0) = TDCR->getBoundaryConditions(stage_num);
        fjac.block<6, 1>(0, i) = MathUtils::forwardFiniteDifferences(f.block<6, 1>(0, 0), f_plus.block<6, 1>(0 , 0), EPS); 

    }

    return 0; 

} 

int LevenbergMarquardtFunctor::values() const 

{

    return m; 

}

int LevenbergMarquardtFunctor::inputs() const 

{

    return n; 

}


LevenbergMarquardtFunctor::~LevenbergMarquardtFunctor() {};