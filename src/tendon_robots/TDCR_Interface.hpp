#include <LevenbergMarquardtFunctor.hpp>
#include <MultistageTDCR_Solver.hpp>
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>

#define assertm(exp, msg) assert(((void)msg, exp))

class TDCR_Interface 

{

    public: 

        TDCR_Interface(MultistageTDCR_Solver& TDCR_);
        void solveForwardKinematics(Eigen::MatrixXd tau, bool print_level);
        std::vector<Eigen::MatrixXd> getJacobians(); 
        virtual ~TDCR_Interface(); 

    private:

        MultistageTDCR_Solver* TDCR;
        LevenbergMarquardtFunctor LMFunctor;

};