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
        void setScaleLoadCell(double value); 
        std::vector<Eigen::MatrixXd> getJacobians(); 
        std::vector<Eigen::MatrixXd> getJacobiansEta(); 
        Eigen::MatrixXd getHighLevelControl(Eigen::Matrix<double, 7, 1> desiredPose);
        Eigen::MatrixXd getHighLevelControl(Eigen::MatrixXd desiredPose);
        Eigen::MatrixXd getDesiredTensions();
        Eigen::MatrixXd getCustomPoseError();
        Eigen::Matrix<double, 7, 1> getPoseError();
        bool checkBoundaryConditions();
        void simulateStep(Eigen::MatrixXd);
        void setWeightsAllStages(Eigen::MatrixXd W_all);
        void setWeightsTip(Eigen::MatrixXd W);
        void trackMeasuredTension(Eigen::MatrixXd tau);
        void setDimensions(double numControlStates, std::vector<Eigen::MatrixXi> CSM);
        virtual ~TDCR_Interface(); 

    private:

        std::vector<Eigen::MatrixXi> controlStatesMatrix;
        Eigen::MatrixXd desiredTensions;
        Eigen::MatrixXd deltaTension;
        MultistageTDCR_Solver* TDCR;
        LevenbergMarquardtFunctor LMFunctor;
        double lambda, Kp, scaleLoadCell;
        Eigen::MatrixXd I_mxm;
        Eigen::Matrix<double, 7, 1> poseTip;
        Eigen::Matrix<double, 7, 1> poseTipError;
        Eigen::MatrixXd weightPriorityTip;
        Eigen::MatrixXd weightPriorityCustom;
        Eigen::MatrixXd customJacobianEta; 
        Eigen::MatrixXd customPose;
        Eigen::MatrixXd JacobianEta; 
        Eigen::MatrixXd customPoseError; 
        bool dimensionsSet; 

};