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

        double simulatedTime;
        TDCR_Interface* tendonRobot;
        Eigen::MatrixXd desiredPose;
        Eigen::MatrixXd controlInput;
        Eigen::MatrixXd desiredTension;
        std_msgs::Float64MultiArray desiredTensionMsg;
        geometry_msgs::PoseStamped poseMsg;

    public: 

        

        tendon_robot_simulator(ros::NodeHandle *nh, TDCR_Interface* tendonRobot_) 
        
        {

            
            boost::thread thread(&tendon_robot_simulator::runDynamicReconfigure, this);
            pub_tensions = nh->advertise<std_msgs::Float64MultiArray>("/desired_tensions", 10);
            pub_pose = nh->advertise<geometry_msgs::PoseStamped>("/tdcr_pose", 10);
            tendonRobot = tendonRobot_;
            desiredPose.resize(tendonRobot->getCustomPose().rows(), tendonRobot->getCustomPose().cols());
            desiredPose << 0.0, 0.0, 1.0, 0.0, 0.0, 0.0;
            simulatedTime = 0.0; 
            simulationLoop();  
            thread.join();         
            ROS_INFO("Closing node...");

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

        void simulationLoop() 
        
        {

            while (ros::ok()) 
            
            {

                controlInput = tendonRobot->getHighLevelControl(desiredPose);
                tendonRobot->simulateStep(controlInput);
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
            desiredPose << msg.px3, msg.py3, eta.w(), eta.vec(); 

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
    stage_tendons << 1, 1, 0, 0, 0, 0,
                     0, 0, 1, 1, 0, 0,
                     0, 0, 0, 0, 1, 1;

    Eigen::MatrixXd routing; 
    routing.resize(3, 6);
    routing.row(2).setZero(); 

    routing(0, 0) = 0.0303;
    routing(1, 0) = 0.0175;
    routing(0, 1) = 0.0303;
    routing(1, 1) = -0.0175;
    routing(0, 2) = -0.0193;
    routing(1, 2) = 0.0230;
    routing(0, 3) = -0.0193;
    routing(1, 3) = -0.0230;
    routing(0, 4) = 0;
    routing(1, 4) = 0.025;
    routing(0, 5) = 0;
    routing(1, 5) = -0.025;


    // Create solver and interfaces.

    MultistageTDCR_Solver tendon_robot(20, 6, 3, i, is, stage_tendons, routing);
    ControllerInterface controller(nlpcapsule);
    TDCR_Interface c(&tendon_robot, &controller); 

    Eigen::MatrixXd tau(6, 1); 
    tau << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; 
    c.solveForwardKinematics(tau, true);
    c.simulateStep(tau);

    Eigen::MatrixXi w1(6, 1);
    Eigen::MatrixXd controlInput; 
    std::vector<Eigen::MatrixXi> CSM; 
    Eigen::MatrixXi CS(1, 1); 
    CS << 2;
    w1 << 0, 1, 3, 4, 5, 6; 
    CSM.push_back(w1); 
    c.setDimensions(6, CSM, CS);

    ros::init(argc, argv, "tdcr_crt");    
    ros::NodeHandle n;
    ros::Rate loop_rate(1000);
    tendon_robot_simulator test(&n, &c);

    free(capsule1);
    free(capsule1step);
    free(capsule2);
    free(capsule2step);

    return 0; 

}