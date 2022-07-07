#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>

int main() 

{

    /*  Convert homogeneous matrices from Optitrack to any specific pose you need. 
    */

    Eigen::Quaterniond q1(2, 0, 1, -3); // Create dummy orientations and normalize them.
    q1.normalize();
    Eigen::Quaterniond q2(1, 0.6, 0, 0); 
    q2.normalize();
    Eigen::Quaterniond q3; 
    Eigen::Quaterniond q_tmp(0, 0, 0, 0);

    Eigen::Vector3d p1(1, 0.5, 0.8);
    Eigen::Vector3d p2(0.7, 0.5, 1.1);
    Eigen::Vector3d p3; 

    q_tmp.vec() = p1 - p2;
    
    q3 = q1.inverse() * q2; // Apply rotation to another rotation.
    q_tmp = q1.inverse() * q_tmp * q1; // Apply rotation to vector. 
    p3 = q_tmp.vec(); 
    
    std::cout << "Check if this is identity: " << q1.toRotationMatrix().inverse() * q2.toRotationMatrix() * q3.toRotationMatrix().inverse() << "\n";
    std::cout << "Check if this adds up to a zero Vector3d: " << q1.toRotationMatrix().inverse() * (p1 - p2) - p3 << "\n";

    return 0;

}