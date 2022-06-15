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
        void setWeights();
        std::vector<Eigen::MatrixXd> getJacobians(); 
        Eigen::MatrixXd getHighLevelControl(Eigen::Matrix<double, 7, 1> desiredPose);
        Eigen::MatrixXd getHighLevelControl(std::vector<Eigen::Matrix<double, 7, 1>> desiredPose);
        void trackMeasuredTension(Eigen::MatrixXd tau);
        virtual ~TDCR_Interface(); 

    private:

        Eigen::MatrixXd desiredTensions;
        Eigen::MatrixXd deltaTension;
        MultistageTDCR_Solver* TDCR;
        LevenbergMarquardtFunctor LMFunctor;

};