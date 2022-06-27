#include "LevenbergMarquardtFunctor.hpp"
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>

LevenbergMarquardtFunctor::LevenbergMarquardtFunctor(MultistageTDCR_Solver& TDCR_)

{

    TDCR = &TDCR_;
    m = TDCR->getNumStages()*6; 
    n = TDCR->getNumStages()*6;

}

LevenbergMarquardtFunctor::LevenbergMarquardtFunctor() 

{}


int LevenbergMarquardtFunctor::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const 

{


    for (unsigned int i = 0; i < TDCR->getNumStages(); ++i) 
    
    {

        TDCR->setInitialConditions(i, x.segment<6>(i*6));
        TDCR->integrateStatesandUpdate(i); 
        fvec.segment<6>(i*6) = TDCR->getBoundaryConditions(i);

    }

    return 0;

}

int LevenbergMarquardtFunctor::df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const 

{

    Eigen::MatrixXd f, f_plus; 
    

    f.resize(TDCR->getNumStages()*6, 1);

    for (unsigned int m = 0; m < TDCR->getNumStages(); ++m) 

    {

        TDCR->setInitialConditions(m, x.segment<6>(m*6));
        TDCR->integrateStatesandUpdate(m); 
        f.block<6, 1>(6*m, 0) = TDCR->getBoundaryConditions(m); 

    }

    f_plus.resize(TDCR->getNumStages()*6, 1);

    for (unsigned int i = 0; i < TDCR->getNumStages()*6; ++i) 

    {

        Eigen::VectorXd x_plus(x);
        x_plus(i) += EPS;

        for (unsigned int k = 0; k < TDCR->getNumStages(); ++k) 
        
        {

            for (unsigned int l = 0; l < TDCR->getNumStages(); ++l) 
            
            {

                TDCR->setInitialConditions(l, x_plus.segment<6>(l*6));
            
            }
            
            TDCR->integrateStatesandUpdate(k);           
            f_plus.block<6, 1>(6*k , 0) = TDCR->getBoundaryConditions(k);
            fjac.block<6, 1>(6*k, i) = MathUtils::forwardFiniteDifferences(f.block<6, 1>(6*k, 0), f_plus.block<6, 1>(6*k , 0), EPS); 

        }

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