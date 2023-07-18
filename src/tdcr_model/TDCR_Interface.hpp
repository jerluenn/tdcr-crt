#include <LevenbergMarquardtFunctor.hpp>
#include <MultistageTDCR_Solver.hpp>
#include <Eigen/Eigen>
#include <unsupported/Eigen/NonLinearOptimization>
#include <ControllerInterface.hpp>

#define assertm(exp, msg) assert(((void)msg, exp))

class TDCR_Interface 

{

    public: 

        TDCR_Interface(MultistageTDCR_Solver* TDCR_, ControllerInterface* MPC_);
        void solveForwardKinematics(Eigen::MatrixXd tau, bool print_level);
        void setScaleLoadCell(double value); 
        std::vector<Eigen::MatrixXd> getJacobians(); 
        std::vector<Eigen::MatrixXd> getJacobiansEta(); 
        Eigen::MatrixXd getHighLevelControl(Eigen::Matrix<double, 7, 1> desiredPose);
        Eigen::MatrixXd getHighLevelControl(Eigen::MatrixXd desiredPose);
        Eigen::MatrixXd getHighLevelControl(Eigen::MatrixXd desiredPose, Eigen::MatrixXd currentPose);
        Eigen::MatrixXd getDesiredTensions();
        Eigen::MatrixXd getCustomPoseError();
        Eigen::MatrixXd getCustomPose();
        Eigen::MatrixXd getTau();
        unsigned int getNumTendons();
        Eigen::MatrixXd saveFullStates(std::string fileName);
        Eigen::MatrixXd getInitialCondition();
        std::vector<Eigen::Matrix<double, 7, 1>> getPoseWorld();
        Eigen::Matrix<double, 7, 1> getPoseError();
        bool checkBoundaryConditions();
        void simulateStep(Eigen::MatrixXd);
        void trackMeasuredTension(Eigen::MatrixXd tau);
        void setDimensions(double numControlStates, std::vector<Eigen::MatrixXi> CSM, Eigen::MatrixXi stagesControlled_);
        virtual ~TDCR_Interface(); 
        MathUtils::Timer timer;

    private:

        std::vector<Eigen::MatrixXi> controlStatesMatrix;
        Eigen::MatrixXi stagesControlled;
        Eigen::MatrixXd desiredTensions;
        Eigen::MatrixXd deltaTension;
        ControllerInterface* MPC;
        MultistageTDCR_Solver* TDCR;
        LevenbergMarquardtFunctor LMFunctor;
        double scaleLoadCell;
        Eigen::MatrixXd I_mxm;
        Eigen::Matrix<double, 7, 1> poseTip;
        Eigen::Matrix<double, 7, 1> poseTipError;
        Eigen::MatrixXd customJacobianEta; 
        Eigen::MatrixXd customPose;
        Eigen::MatrixXd JacobianEta; 
        Eigen::MatrixXd customPoseError; 

        bool dimensionsSet; 

};