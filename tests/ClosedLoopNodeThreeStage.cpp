#include <TDCR_Interface.hpp>
#include "acados_sim_solver_multistage_straight_integrator1.h"
#include "acados_sim_solver_multistage_straight_integrator2.h"
#include "acados_sim_solver_multistage_straight_integrator3.h"
#include "acados_sim_solver_multistage_straight_step_integrator1.h"
#include "acados_sim_solver_multistage_straight_step_integrator2.h"
#include "acados_sim_solver_multistage_straight_step_integrator3.h"
#include "acados_solver_tdcr_lmpc.h"
#include "ros/ros.h"
#include <boost/thread.hpp>
#include "boost/thread/mutex.hpp"
#include <tdcr_crt/SimulationParamsConfig.h>
#include "dynamic_reconfigure/server.h"
#include "std_msgs/Float64MultiArray.h"
#include "std_msgs/Bool.h"
#include "geometry_msgs/PoseStamped.h"

#define PI 3.14159 

class tendon_robot_simulator 

{

    private: 

        ros::Publisher pub_tensions;
        ros::Publisher pub_pose; 
        ros::Subscriber sub_current_pose_03; 
        ros::Subscriber sub_current_pose_02; 
        ros::Subscriber sub_current_pose_01; 
        ros::Subscriber sub_continuum_robot_tension; 

        double simulatedTime;
        TDCR_Interface* tendonRobot;
        Eigen::MatrixXd desiredPose;
        Eigen::MatrixXd controlInput;
        Eigen::MatrixXd tensionInput_updateJacobians;
        Eigen::MatrixXd desiredTension;
        Eigen::Matrix<double, 7, 1> pose01;
        Eigen::Matrix<double, 7, 1> pose02;
        Eigen::Matrix<double, 7, 1> pose03; 
        Eigen::MatrixXd currentPose;
        std_msgs::Float64MultiArray desiredTensionMsg;
        geometry_msgs::PoseStamped poseMsg;

    public: 

        

        tendon_robot_simulator(ros::NodeHandle *nh, TDCR_Interface* tendonRobot_) 
        
        {

            boost::thread thread(&tendon_robot_simulator::runDynamicReconfigure, this);
            pub_tensions = nh->advertise<std_msgs::Float64MultiArray>("/desired_tensions", 10);
            pub_pose = nh->advertise<geometry_msgs::PoseStamped>("/tdcr_pose", 10);
            sub_continuum_robot_tension = nh->subscribe("/distal_tension", 1, &tendon_robot_simulator::continuum_robot_tension_callback , this);
            sub_current_pose_03 = nh->subscribe("/relative_pose_03", 1, &tendon_robot_simulator::relativePose03_callback , this);
            sub_current_pose_02 = nh->subscribe("/relative_pose_02", 1, &tendon_robot_simulator::relativePose02_callback , this);
            sub_current_pose_01 = nh->subscribe("/relative_pose_01", 1, &tendon_robot_simulator::relativePose01_callback , this);
            tendonRobot = tendonRobot_;
            pose01.setZero();
            pose02.setZero();
            pose03.setZero();
            tensionInput_updateJacobians.resize(tendonRobot->getNumTendons(), 1);
            tensionInput_updateJacobians.setZero();
            desiredPose.resize(tendonRobot->getCustomPose().rows(), tendonRobot->getCustomPose().cols()); 
            // desiredPose << 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
            currentPose.resize(tendonRobot->getCustomPose().rows() + tendonRobot->getNumTendons(), tendonRobot->getCustomPose().cols());
            currentPose.setZero(); 
            desiredPose << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
            simulatedTime = 0.0; 
            controlLoop();  
            thread.join();         
            ROS_INFO("Closing node...");

        }

        void continuum_robot_tension_callback(const std_msgs::Float64MultiArray& msg) {

            tensionInput_updateJacobians = Eigen::VectorXd::Map(msg.data.data(), msg.data.size());

        } 

        void relativePose02_callback(const geometry_msgs::PoseStamped::ConstPtr& relative_pose) 
        
        {

            pose02(0, 0) = relative_pose->pose.position.x;
            pose02(1, 0) = relative_pose->pose.position.y;
            pose02(2, 0) = relative_pose->pose.position.z;
            pose02(3, 0) = relative_pose->pose.orientation.w;
            pose02(4, 0) = relative_pose->pose.orientation.x;
            pose02(5, 0) = relative_pose->pose.orientation.y;
            pose02(6, 0) = relative_pose->pose.orientation.z;


        } 

        void relativePose01_callback(const geometry_msgs::PoseStamped::ConstPtr& relative_pose) 
        
        {

            pose01(0, 0) = relative_pose->pose.position.x;
            pose01(1, 0) = relative_pose->pose.position.y;
            pose01(2, 0) = relative_pose->pose.position.z;
            pose01(3, 0) = relative_pose->pose.orientation.w;
            pose01(4, 0) = relative_pose->pose.orientation.x;
            pose01(5, 0) = relative_pose->pose.orientation.y;
            pose01(6, 0) = relative_pose->pose.orientation.z;


        } 

        void relativePose03_callback(const geometry_msgs::PoseStamped::ConstPtr& relative_pose) 
        
        {

            pose03(0, 0) = relative_pose->pose.position.x;
            pose03(1, 0) = relative_pose->pose.position.y;
            pose03(2, 0) = relative_pose->pose.position.z;
            pose03(3, 0) = relative_pose->pose.orientation.w;
            pose03(4, 0) = relative_pose->pose.orientation.x;
            pose03(5, 0) = relative_pose->pose.orientation.y;
            pose03(6, 0) = relative_pose->pose.orientation.z;


        } 

        void publishDesiredTension() 
        
        {

            desiredTension = tendonRobot->getDesiredTensions();

            std::vector<double> desiredTensionVec(desiredTension.data(), desiredTension.data() + desiredTension.size());
            desiredTensionMsg.data.clear();
            desiredTensionMsg.data.insert(desiredTensionMsg.data.end(), desiredTensionVec.begin(), desiredTensionVec.end());
            pub_tensions.publish(desiredTensionMsg);

        }

        void publishPose() 
        
        {

            poseMsg.pose.position.x = tendonRobot->getPoseWorld()[2](0, 0);
            poseMsg.pose.position.y = tendonRobot->getPoseWorld()[2](1, 0);
            poseMsg.pose.position.z = tendonRobot->getPoseWorld()[2](2, 0);
            poseMsg.pose.orientation.w = tendonRobot->getPoseWorld()[2](3, 0);  
            poseMsg.pose.orientation.x = tendonRobot->getPoseWorld()[2](4, 0);
            poseMsg.pose.orientation.y = tendonRobot->getPoseWorld()[2](5, 0);
            poseMsg.pose.orientation.z = tendonRobot->getPoseWorld()[2](6, 0);

            pub_pose.publish(poseMsg);

        }

        std::string toString(const Eigen::MatrixXd& mat) 
        
        {

            std::stringstream ss; 
            ss << mat; 
            return ss.str();

        }

       void updateCurrentPose() 
       
       {

            // int numStates = 6; 

            // currentPose(0, 0) = pose03(0, 0); 
            // currentPose(1, 0) = pose03(1, 0);
            // currentPose(2, 0) = pose03(3, 0);
            // currentPose(3, 0) = pose03(4, 0);
            // currentPose(4, 0) = pose03(5, 0);
            // currentPose(5, 0) = pose03(6, 0);

            // for (unsigned int i = 0; i < tendonRobot->getNumTendons(); i++) 
            
            // {

            //     if (tendonRobot->getTau()(i, 0) < 0) 
                
            //     {
            //         currentPose(numStates + i) = 0.01; 

            //     }


            //     else 
                
            //     {

            //         currentPose(numStates + i) = tendonRobot->getTau()(i, 0);

            //     }

            // }


            int numStates = 10; 

            currentPose(0, 0) = pose01(0, 0); 
            currentPose(1, 0) = pose01(1, 0);
            currentPose(2, 0) = pose02(0, 0);
            currentPose(3, 0) = pose02(1, 0);
            currentPose(4, 0) = pose03(0, 0);
            currentPose(5, 0) = pose03(1, 0);
            currentPose(6, 0) = pose03(3, 0);
            currentPose(7, 0) = pose03(4, 0);
            currentPose(8, 0) = pose03(5, 0);
            currentPose(9, 0) = pose03(6, 0);

            for (unsigned int i = 0; i < tendonRobot->getNumTendons(); i++) 
            
            {

                if (tendonRobot->getTau()(i, 0) < 0) 
                
                {
                    currentPose(numStates + i) = 0.01; 

                }


                else 
                
                {

                    currentPose(numStates + i) = tendonRobot->getTau()(i, 0);

                }

            }

       }
       
        void controlLoop() 
        
        {

            while (ros::ok()) 
            
            {

                tendonRobot->timer.tic();
                updateCurrentPose();
                tendonRobot->trackMeasuredTension(tensionInput_updateJacobians);
                controlInput = tendonRobot->getHighLevelControl(desiredPose, currentPose);
                tendonRobot->timer.toc();
                tendonRobot->checkBoundaryConditions();
                publishPose();
                publishDesiredTension();
                ros::spinOnce();
               
            }


        }

        void runDynamicReconfigure() 

        {

		ROS_INFO("Setting up the dynamic reconfigure panel and server");

			dynamic_reconfigure::Server<tdcr_crt::SimulationParamsConfig> server;
			dynamic_reconfigure::Server<tdcr_crt::SimulationParamsConfig>::CallbackType f;
			f = boost::bind(&tendon_robot_simulator::callbackDynamicReconfigure, this, _1, _2);
			server.setCallback(f);

		    ros::spin();        

        }

        void callbackDynamicReconfigure(tdcr_crt::SimulationParamsConfig& msg, uint32_t level) 

        {

            Eigen::Quaterniond eta; 
            Eigen::Matrix<double, 3, 1> eulerAngles(msg.phi, msg.theta, msg.psi);
            eta = MathUtils::eul2quat_deg(eulerAngles);
            desiredPose << msg.px1, msg.py1, msg.px2, msg.py2, msg.px3, msg.py3, eta.w(), eta.vec();
            // desiredPose << msg.px3, msg.py3, eta.w(), eta.vec(); 

        }


};

int main(int argc, char **argv) 

{

    // Prepare acados integrators here.

    sim_solver_capsule *capsule1 = multistage_straight_integrator1_acados_sim_solver_create_capsule();
    multistage_straight_integrator1_acados_sim_create(capsule1);
    sim_solver_capsule *capsule2 = multistage_straight_integrator2_acados_sim_solver_create_capsule();
    multistage_straight_integrator2_acados_sim_create(capsule2);
    sim_solver_capsule *capsule3 = multistage_straight_integrator2_acados_sim_solver_create_capsule();
    multistage_straight_integrator3_acados_sim_create(capsule3);
    sim_solver_capsule *capsule1step = multistage_straight_integrator1_acados_sim_solver_create_capsule();
    multistage_straight_step_integrator1_acados_sim_create(capsule1step);
    sim_solver_capsule *capsule2step = multistage_straight_integrator2_acados_sim_solver_create_capsule();
    multistage_straight_step_integrator2_acados_sim_create(capsule2step);
    sim_solver_capsule *capsule3step = multistage_straight_integrator2_acados_sim_solver_create_capsule();
    multistage_straight_step_integrator3_acados_sim_create(capsule3step);
    tdcr_lmpc_solver_capsule *nlpcapsule = tdcr_lmpc_acados_create_capsule();
    tdcr_lmpc_acados_create(nlpcapsule);

    std::vector<IntegrationInterface> i;
    std::vector<IntegrationInterface> is;

    IntegrationInterface i1(capsule1), i2(capsule2), ii1(capsule1step), ii2(capsule2step), i3(capsule3), ii3(capsule3step); 
    i.push_back(i1);
    i.push_back(i2);
    i.push_back(i3);
    is.push_back(ii1);
    is.push_back(ii2);
    is.push_back(ii3);

    // Set routing here.

    Eigen::MatrixXd stage_tendons;
    stage_tendons.resize(3, 6);
    stage_tendons << 0, 0, 1, 0, 0, 1,
                     0, 1, 0, 0, 1, 0,
                     1, 0, 0, 1, 0, 0;

    Eigen::MatrixXd routing; 
    routing.resize(3, 6);
    routing.row(2).setZero(); 

    routing(0, 0) = 0.0189;
    routing(1, 0) = 0.012164;
    routing(0, 1) = -0.003202;
    routing(1, 1) = 0.0223;
    routing(0, 2) = -0.0189;
    routing(1, 2) = 0.012164;
    routing(0, 3) = 0.0189;
    routing(1, 3) = -0.012164;
    routing(0, 4) = -0.003202;
    routing(1, 4) = -0.0223;
    routing(0, 5) = -0.0189;
    routing(1, 5) = -0.012164;

    // [[0.01351, 0.0210331, 0], [0.0035578, 0.024745, 0], [-0.016375, 0.0188937, 0], [0.01351, -0.0210331, 0], [0.0035578, -0.024745, 0], [-0.016375, -0.0188937, 0]]

    // Create solver and interfaces.

    MultistageTDCR_Solver tendon_robot(20, 6, 3, i, is, stage_tendons, routing);
    ControllerInterface controller(nlpcapsule);
    TDCR_Interface c(&tendon_robot, &controller); 

    Eigen::MatrixXd tau(6, 1); 
    tau << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; 
    c.solveForwardKinematics(tau, true);
    c.simulateStep(tau);

    Eigen::MatrixXi w1(2, 1), w2(2, 1), w3(6, 1);
    Eigen::MatrixXd controlInput; 
    std::vector<Eigen::MatrixXi> CSM; 
    Eigen::MatrixXi CS(3, 1); 
    CS << 0, 1, 2;
    w1 << 0, 1;
    w2 << 0, 1;
    w3 << 0, 1, 3, 4, 5, 6; 
    CSM.push_back(w1); 
    CSM.push_back(w2); 
    CSM.push_back(w3); 
    c.setDimensions(10, CSM, CS);

    // Eigen::MatrixXi w3(6, 1);
    // Eigen::MatrixXd controlInput; 
    // std::vector<Eigen::MatrixXi> CSM; 
    // Eigen::MatrixXi CS(1, 1); 
    // CS << 2;
    // w3 << 0, 1, 3, 4, 5, 6; 
    // CSM.push_back(w3); 
    // c.setDimensions(6, CSM, CS);

    ros::init(argc, argv, "tdcr_crt_closed_loop");    
    ros::NodeHandle n;
    ros::Rate loop_rate(1000);
    tendon_robot_simulator test(&n, &c);

    free(capsule1);
    free(capsule1step);
    free(capsule2);
    free(capsule2step);

    return 0; 

}