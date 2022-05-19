#ifndef LevenbergMarquardtFunctor_H
#define LevenbergMarquardtFunctor_H

#include <MultistageTDCR_Solver.hpp>

class LevenbergMarquardtFunctor

{

    public: 

        LevenbergMarquardtFunctor(MultistageTDCR_Solver& TDCR_);
        LevenbergMarquardtFunctor();
        virtual ~LevenbergMarquardtFunctor();
        void solveForwardKinematics();

        // Compute 'm' errors, one for each data point, for the given parameter values in 'x'
        int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const;

        // Compute the jacobian of the errors
        int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const;

        // Number of data points, i.e. values.
        int m;

        // Returns 'm', the number of values.
        int values() const ;

        // The number of parameters, i.e. inputs.
        int n;

        // Returns 'n', the number of inputs.
        int inputs() const ;

    private: 

        MultistageTDCR_Solver* TDCR; 
        static constexpr double EPS = 1e-9;

};

#endif