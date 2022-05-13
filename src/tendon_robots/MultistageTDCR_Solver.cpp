#include "MultistageTDCR_Solver.hpp"

Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");

MultistageTDCR_Solver::MultistageTDCR_Solver(int numTendons, int numStages, std::vector<IntegrationInterface> integrators_, std::vector<IntegrationInterface> integratorsStep_, Eigen::MatrixXd stage_tendons, Eigen::MatrixXd routing_)
{               

    v << 0.0, 0.0, 1.0;
    integrators = integrators_;
    integratorsStep = integratorsStep_;
    setNumTendons(numTendons, routing_);
    setNumStages(numStages);
    initialiseJacobianMatrices(stage_tendons);

    testFunction();


}

void MultistageTDCR_Solver::convertStageTendonsIndex() 

{

    std::vector<int> indices; 
    Eigen::VectorXi indices_eig; 
    int vector_size = 0; 
    int index = 0; 

    for (auto it : stageTendons.rowwise()) 

        {

        for (int i = 0; i < it.size();  ++i) {

            if (it(i) == 1) {

                indices.push_back(i);

            } 

        
        }

        indices_eig = Eigen::VectorXi::Map(indices.data(), indices.size());
        stageTendonsIndex_Termination.push_back(indices_eig);
        indices.clear();
        
        }


    for (unsigned int j = 0; j < num_stages; ++j) 

    {
        
        for (int k = num_stages - j - 1; k >= 0; --k)

        {

            vector_size += stageTendonsIndex_Termination[k].size();

        }

        indices_eig.resize(vector_size);

        for (int l = j; l < int(num_stages); ++l)

        {


            indices_eig.segment(index, stageTendonsIndex_Termination[l].size()) = stageTendonsIndex_Termination[l];
            
            
            index += stageTendonsIndex_Termination[l].size();

        }

        stageTendonsIndex_Jacobians.push_back(indices_eig);
        vector_size = 0; 
        index = 0;

    } 

}

void MultistageTDCR_Solver::initialiseJacobianMatrices(Eigen::MatrixXd stage_tendons) {

    assertm(stage_tendons.rows() == num_stages, "Number of rows in stage tendons must equal number of stages!");
    assertm(stage_tendons.cols() == num_tendons, "Number of cols in stage tendons must equal number of tendons!");
    stageTendons = stage_tendons;    
    Eigen::MatrixXd tmp; 

    convertStageTendonsIndex();
    AdjointMatrix.setIdentity();

    for(unsigned int i = 0; i < num_stages; ++i) {

        tmp.resize(6,6); 
        tmp.setZero();
        B_yu_Nplus1.push_back(tmp);
        tmp.resize(6,num_tendons);
        tmp.setZero();
        J_q.push_back(tmp);
        tmp.resize(6,num_tendons);
        tmp.setZero();
        J_b.push_back(tmp);
        tmp.resize(6,num_tendons);
        tmp.setZero();
        E_q.push_back(tmp);
        tmp.resize(6,num_tendons);
        tmp.setZero();
        B_q.push_back(tmp);
        tmp.setZero();
        tmp.resize(6,6);
        B_yu.push_back(tmp);
        tmp.resize(6,6);
        tmp.setZero();
        E_yu.push_back(tmp);
        tmp.resize(6,num_tendons);
        tmp.setZero();
        J_world.push_back(tmp);
        

    }

}


Eigen::MatrixXd MultistageTDCR_Solver::integrateWithIncrement(unsigned int index, unsigned int stage_num) {

    x = initialConditions[stage_num];
    x(index) += EPS;
    return integrators[stage_num].integrate(x);

}

Eigen::MatrixXd MultistageTDCR_Solver::addIncrement(unsigned int index, unsigned int stage_num) 

{

    x = initialConditions[stage_num]; 
    x(index) += EPS; 
    return x; 

}

Eigen::Matrix<double, 6, 1> MultistageTDCR_Solver::getBoundaryConditions(unsigned int stage_num, Eigen::MatrixXd robotStates_, Eigen::MatrixXd initialConditions_) 

{

    // Initial Conditions here must be at stage N + 1!

    Eigen::Matrix<double, 6, 1> boundaryConditions;
    Eigen::Matrix<double, 6, 1> pointForceMoment;
    Eigen::Matrix<double, 6, 1> internalForcesandMoments_N_l; // internal wrench at stage N, s = 0.
    Eigen::Matrix<double, 6, 1> internalForcesandMoments_Nplus1_0; // internal wrench at stage N-1, s = l. 
    Eigen::Matrix<double, 6, 6> R_diag; 
    Eigen::Matrix<double, 3, 3> zeros_3x3; 
    zeros_3x3.setZero();
    R.resize(9,1);
    pointForceMoment.setZero();

    R = robotStates_.block<9,1>(3,0); 
    R.resize(3,3); 
    internalForcesandMoments_N_l = robotStates_.block<6, 1>(num_p + num_R, 0);

    if (stage_num == num_stages - 1) 
    {

        internalForcesandMoments_Nplus1_0.setZero(); 

    }

    else 

    {

        internalForcesandMoments_Nplus1_0 = initialConditions_.block<6, 1>(num_p + num_R, 0);

    }

    pointForceMoment = getPointForceMoment(stage_num);
    R_diag << R, zeros_3x3, zeros_3x3, R; 

    // boundaryConditions = internalForcesandMoments_N_l - R_diag * internalForcesandMoments_Nplus1_0 - pointForceMoment;
    boundaryConditions = - internalForcesandMoments_N_l + R_diag * internalForcesandMoments_Nplus1_0 + pointForceMoment;

    return boundaryConditions;

} 

Eigen::Matrix<double, 6, 1> MultistageTDCR_Solver::getBoundaryConditions(unsigned int stage_num) 

{

    // Initial Conditions here must be at stage N + 1!

    Eigen::Matrix<double, 6, 1> boundaryConditions;
    Eigen::Matrix<double, 6, 1> pointForceMoment;
    Eigen::Matrix<double, 6, 1> internalForcesandMoments_N_l; // internal wrench at stage N, s = 0.
    Eigen::Matrix<double, 6, 1> internalForcesandMoments_Nplus1_0; // internal wrench at stage N-1, s = l. 
    Eigen::Matrix<double, 6, 6> R_diag; 
    Eigen::Matrix<double, 3, 3> zeros_3x3; 
    zeros_3x3.setZero();
    R.resize(9,1);
    pointForceMoment.setZero();
    R = robotStates[stage_num].block<9,1>(3,0); 
    R.resize(3,3);

    internalForcesandMoments_N_l = robotStates[stage_num].block<6, 1>(num_p + num_R, 0);

    if (stage_num == num_stages - 1) 
    {

        internalForcesandMoments_Nplus1_0.setZero(); 

    }

    else 

    {

        internalForcesandMoments_Nplus1_0 = initialConditions[stage_num + 1].block<6, 1>(num_p + num_R, 0);

    }

    pointForceMoment = getPointForceMoment(stage_num);

    R_diag << R, zeros_3x3, zeros_3x3, R; 


    // boundaryConditions = internalForcesandMoments_N_l - R_diag * internalForcesandMoments_Nplus1_0 - pointForceMoment;
    boundaryConditions = - internalForcesandMoments_N_l + R_diag * internalForcesandMoments_Nplus1_0 + pointForceMoment;

    return boundaryConditions;

} 

Eigen::Matrix<double, 6, 1> MultistageTDCR_Solver::getBoundaryConditions(unsigned int stage_num, Eigen::MatrixXd robotStates_) 

{

    // Initial Conditions here must be at stage N + 1!

    Eigen::Matrix<double, 6, 1> boundaryConditions;
    Eigen::Matrix<double, 6, 1> pointForceMoment;
    Eigen::Matrix<double, 6, 1> internalForcesandMoments_N_l; // internal wrench at stage N-1, s = l.
    Eigen::Matrix<double, 6, 1> internalForcesandMoments_Nplus1_0;  // internal wrench at stage N, s = 0.
    Eigen::Matrix<double, 6, 6> R_diag; 
    Eigen::Matrix<double, 3, 3> zeros_3x3; 
    zeros_3x3.setZero();
    pointForceMoment.setZero();

    internalForcesandMoments_N_l = robotStates_.block<6, 1>(num_p + num_R, 0);

    if (stage_num == num_stages - 1) 
    {

        internalForcesandMoments_Nplus1_0.setZero(); 

    }

    else 

    {

        internalForcesandMoments_Nplus1_0 = initialConditions[stage_num + 1].block<6, 1>(num_p + num_R, 0);

    }

    pointForceMoment = getPointForceMoment(stage_num, robotStates_);

    R.resize(9,1);
    
    R = robotStates_.block<9,1>(3,0); 
    R.resize(3,3); 

    R_diag << R, zeros_3x3, zeros_3x3, R;

    boundaryConditions = - internalForcesandMoments_N_l + R_diag * internalForcesandMoments_Nplus1_0 + pointForceMoment;

    return boundaryConditions;

} 

Eigen::Matrix<double, 6, 1> MultistageTDCR_Solver::getPointForceMoment(unsigned int stage_num) 

{

    /* Get Point Force and Point Moment in the local frame. (Frame N) */

    Eigen::Matrix<double, 6, 1> pointForceMoment;
    R.resize(9,1);
    pointForceMoment.setZero();
    Eigen::Vector3d PointForce;
    Eigen::Vector3d tmpVector; 


    for (auto k : stageTendonsIndex_Termination[stage_num]) 
    {

        R = robotStates[stage_num].block<9,1>(3,0); 
        R.resize(3,3);
        PointForce = R*v*initialConditions[stage_num].block(num_p + num_R + num_m + num_n + k, 0, 1, 1);
        pointForceMoment.block<3, 1>(0, 0) -= PointForce;// Point Force.
        tmpVector = (R*routing.col(k));
        pointForceMoment.block<3, 1>(3, 0) -= tmpVector.cross(PointForce); // Point Moment
        
    }
              
    return pointForceMoment;

} 


Eigen::Matrix<double, 6, 1> MultistageTDCR_Solver::getPointForceMoment(unsigned int stage_num, Eigen::MatrixXd robotStates_) 

{

    /* Get Point Force and Point Moment in the local frame. (Frame N) */

    Eigen::Matrix<double, 6, 1> pointForceMoment;
    R.resize(9,1);
    pointForceMoment.setZero();
    Eigen::Vector3d PointForce;
    Eigen::Vector3d tmpVector; 


    for (auto k : stageTendonsIndex_Termination[stage_num]) 
    {

        R = robotStates_.block<9,1>(3,0); 
        R.resize(3,3);
        PointForce = R*v*robotStates_.block(num_p + num_R + num_m + num_n + k, 0, 1, 1);
        pointForceMoment.block<3, 1>(0, 0) -= PointForce;// Point Force.
        tmpVector = (R*routing.col(k));
        pointForceMoment.block<3, 1>(3, 0) -= tmpVector.cross(PointForce); // Point Moment
        
    }
              
    return pointForceMoment;

} 

void MultistageTDCR_Solver::setTau(Eigen::MatrixXd tau) 

{

    assertm(num_tendons == tau.rows(), "tau must be the same size as the number of stages!");

    for (unsigned int i = 0; i < num_stages; ++i) {

        initialConditions[i].block(num_p + num_R + num_m + num_n, 0, num_tendons, 1) = tau;

    }

}


void MultistageTDCR_Solver::setNumTendons(int num, Eigen::MatrixXd routing_) 

{

    num_tendons = num;
    num_total = num_p + num_R + num_m + num_n + num_alpha + num_tendons;

    x.resize(num_total, 1);
    routing.resize(3, num_total);
    routing = routing_;

}

void MultistageTDCR_Solver::setNumStages(int num) 

{

    num_stages = num;
    Eigen::MatrixXd states;
    states.resize(num_total, 1);
    states.setZero();
    states.block<9, 1>(3, 0) << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;

    for (unsigned int i = 0; i < num_stages; ++i ) {

        initialConditions.push_back(states);
        robotStates.push_back(states);

    } 


    assertm(num_stages == integrators.size(), "Checking integrator list size."); 
    assertm(num_stages == integratorsStep.size(), "Checking integrator list size.");

}


std::vector<Eigen::MatrixXd> MultistageTDCR_Solver::getRobotStates(bool print_level, bool csv_level)
{


    if (print_level == true) {

        for (unsigned int i = 0; i < num_stages; ++i) {

            std::cout << "Initial conditions at stage " << i << ": \n";
            std::cout << initialConditions[i].transpose().format(OctaveFmt) << "\n";
            std::cout << "Robot states at stage " << i << ": \n";
            std::cout << robotStates[i].transpose().format(OctaveFmt) << "\n";

        }

    } 

    return robotStates;

}

Eigen::MatrixXd MultistageTDCR_Solver::integrateStates(unsigned int stage_num)
{

    return integrators[stage_num].integrate(initialConditions[stage_num]);

} 

Eigen::MatrixXd MultistageTDCR_Solver::integrateStep(unsigned int stage_num)

{

    return integrators[stage_num].integrate(initialConditions[stage_num]);

}

void MultistageTDCR_Solver::solveJacobians() 

{


    std::vector<Eigen::Matrix<double, 3, 1>> positions; 
    std::vector<Eigen::Matrix3d> rotations; 
    std::vector<Eigen::Matrix<double, 6, 1>> boundaryConditions;

    Eigen::Matrix<double, 3, 1> position;
    
    Eigen::Matrix3d rotation; 
    Eigen::Matrix<double, 6, 1> BC; 

    Eigen::MatrixXd linear_twist; 
    linear_twist.resize(3, 1); 
    Eigen::MatrixXd rotation_twist; 
    rotation_twist.resize(3, 3); 
    Eigen::MatrixXd boundary_twist; 
    boundary_twist.resize(6, 1); 

    Eigen::MatrixXd initCon; 
    initCon.resize(num_total, 1); 

    Eigen::Matrix4d se3; 

    R.resize(9, 1);


    Eigen::MatrixXd robotStates_Increment;
    robotStates_Increment.resize(num_total, 1);

    for (unsigned int i = 0; i < num_stages; ++i) 
    
    {

        position = robotStates[i].block<3, 1>(0,0);
        R = robotStates[i].block<9, 1>(3, 0);
        R.resize(3, 3);

        rotation << R;

        BC = getBoundaryConditions(i, robotStates[i]); 

        positions.push_back(position);
        rotations.push_back(rotation);
        boundaryConditions.push_back(BC);

    } 

    for (int ii = 0; ii < 6; ++ii) 
    
    {

        for (unsigned int k = 0; k < num_stages; ++k) 
        {

            robotStates_Increment = integrateWithIncrement(ii + num_p + num_R, k);
            position = robotStates_Increment.block<3, 1>(0, 0); 
            R = robotStates_Increment.block<9, 1>(3, 0);
            R.resize(3, 3);
            rotation << R; 
            BC = getBoundaryConditions(k, robotStates_Increment);

            linear_twist = MathUtils::forwardFiniteDifferences(positions[k], position, EPS);
            rotation_twist = MathUtils::forwardFiniteDifferences(rotations[k], rotation, EPS);
            
            boundary_twist = MathUtils::forwardFiniteDifferences(boundaryConditions[k], BC, EPS);

            se3.block<3, 3>(0, 0) = rotation_twist*rotations[k].transpose();
            se3.block<3, 1>(0, 3) = linear_twist;

            E_yu[k].col(ii) = MathUtils::se3toVec(se3);
            B_yu[k].col(ii) = boundary_twist; 

        } 
    

        for (unsigned int k = 0; k < num_stages; ++k) 

        {

            if (k + 1 != num_stages) 

            {

                initCon = addIncrement(num_p + num_R + ii, k+1);
                BC = getBoundaryConditions(k, robotStates[k], initCon);
                boundary_twist = MathUtils::forwardFiniteDifferences(boundaryConditions[k], BC, EPS);
                B_yu_Nplus1[k].col(ii) = boundary_twist;

            }

        }

    }

    for (unsigned int j = 0; j < num_stages; ++j) 
    
    {

        for (auto it : stageTendonsIndex_Jacobians[j])

            {

                robotStates_Increment = integrateWithIncrement(num_p + num_R + num_m + num_n + it, j);
                position = robotStates_Increment.block<3, 1>(0, 0); 
                R = robotStates_Increment.block<9, 1>(3, 0);
                R.resize(3, 3);
                rotation << R; 
                BC = getBoundaryConditions(j, robotStates_Increment);              

                linear_twist = MathUtils::forwardFiniteDifferences(positions[j], position, EPS);
                rotation_twist = MathUtils::forwardFiniteDifferences(rotations[j], rotation, EPS);
                boundary_twist = MathUtils::forwardFiniteDifferences(boundaryConditions[j], BC, EPS);

                se3.block<3, 3>(0, 0) = rotation_twist*rotations[j].transpose();
                se3.block<3, 1>(0, 3) = linear_twist;

                E_q[j].col(it) = MathUtils::se3toVec(se3);
                B_q[j].col(it) = boundary_twist;

            }

    } 

    for (int l = num_stages - 1; l >= 0; --l) 

    {

        if (l == int(num_stages - 1))      
        {

            J_b[l] = -B_yu[l].completeOrthogonalDecomposition().pseudoInverse()*B_q[l];
            J_q[l] = E_q[l] + E_yu[l]*J_b[l];

            // std::cout << "Debugging info at stage " << l << ": \n";
            // std::cout << "J_b[l]: " << J_b[l].format(OctaveFmt) << "\n";
            // std::cout << "J_q[l]: " << J_q[l].format(OctaveFmt) << "\n\n\n";

        }

        else 
        {

            J_b[l] = -B_yu[l].completeOrthogonalDecomposition().pseudoInverse()*(B_q[l] + B_yu_Nplus1[l]*J_b[l + 1]);
            J_q[l] = E_q[l] + E_yu[l]*J_b[l];

            // std::cout << "Debugging info at stage " << l << ": \n";
            // std::cout << "B_yu[l+1]" << B_yu[l+1].format(OctaveFmt) << "\n";
            // std::cout << "B_yu_Nplus1" << B_yu_Nplus1[l].format(OctaveFmt) << "\n";
            // std::cout << "R " << R << "\n";
            // std::cout << "B_q[l]" << B_q[l].format(OctaveFmt) << "\n";
            // std::cout << "J_b[l+1]: " << J_b[l+1].format(OctaveFmt) << "\n";
            // std::cout << "J_q[l]: " << J_q[l].format(OctaveFmt) << "\n\n\n";

        }
        
    }

    Eigen::Matrix3d R_tmp; 
    Eigen::Matrix<double, 6, 6> R_6x6;
    R_6x6.setZero();
    R_tmp.setIdentity();
    Eigen::Matrix<double, 3, 1> pos_tmp;
    Eigen::MatrixXd J_world_tmp;
    J_world_tmp.resize(6, num_tendons);
    J_world_tmp.setZero();
    Eigen::MatrixXd rhs_matrix; 
    rhs_matrix.resize(6 , num_tendons);
    rhs_matrix.setZero();

    // Get first J_world_tmp here. 

    J_world_tmp = J_q[0];
    J_world[0] = J_world_tmp;

    for (unsigned int k = 1; k < num_stages; ++k)

    {

        R.resize(9, 1); 
        R = robotStates[k-1].block<9, 1>(3, 0); 
        R.resize(3, 3);
        R_tmp *= R; 
        R_6x6.block<3, 3>(0, 0) = R_tmp; 
        R_6x6.block<3, 3>(3, 3) = R_tmp;
        pos_tmp = robotStates[k].block<3, 1>(0, 0);  
        rhs_matrix.block<3, 6>(0, 0) = -MathUtils::skew_m(R_tmp*pos_tmp) * J_world[k - 1].block(3, 0, 3, num_tendons);

        std::cout << R_tmp << "\n";

        J_world_tmp += R_6x6 * J_q[k] + rhs_matrix;
        J_world[k] = J_world_tmp;

    }  


}

void MultistageTDCR_Solver::solveForwardKinematics()
{

}


void MultistageTDCR_Solver::simulateStep(Eigen::MatrixXd tau)
{

    assertm(tau.cols() == 1, "tau must have only one column, i.e. a vector."); 
    assertm(tau.rows() == B_q[0].cols(), "tau.rows() must match B_q.cols().");

    Eigen::MatrixXd db(6, 1); 

    solveJacobians();

    for (unsigned int i = 0; i < num_stages; ++i) 
    
    {


        db = J_b[i]*tau;
        initialConditions[i].block<6, 1>(num_p + num_R, 0) += db*dt; 
        initialConditions[i].block(num_p + num_R + num_m + num_n, 0, num_tendons, 1) += tau*dt;
        robotStates[i] = integrateStates(i);


    }

    

}

void MultistageTDCR_Solver::testFunction() 

{

    Eigen::Matrix<double, 6, 1> ic; 
    ic << 8.69005567e-01, -1.68328807e-27, -5.43264986e+00, -2.19026323e-27,
        2.73144032e-01, -9.56647789e-29;
    setInitialConditions(0, ic);
    Eigen::Matrix<double, 6, 1> ic1;
    ic1 <<  5.37588000e-01,  4.28851040e-28, -4.88462763e-27, -2.22218785e-27,
        5.36810924e-02, -1.02139239e-28;
    setInitialConditions(1, ic1);
    Eigen::MatrixXd tau;
    tau.resize(6, 1); 
    tau << 5, 0, 0, 0, 0, 0;
    setTau(tau);

    robotStates[0] = integrateStates(0);
    robotStates[1] = integrateStates(1);

    Eigen::MatrixXd tau1; 
    tau1.resize(6, 1);
    tau1 << 0, 1, 1, 2, 0, 0;
    std::vector<Eigen::Matrix<double, 3, 1>> p; 
    Eigen::Matrix<double, 3, 1> p_;
    p_.setZero();
    p.push_back(p_);
    p.push_back(p_);
    std::cout << robotStates[0].block<3, 1>(0, 0) << "\n";
    p[0] = robotStates[0].block<3, 1>(0, 0); 
    p[1] = robotStates[1].block<3, 1>(0, 0);

    Eigen::Matrix<double, 3, 1> p_world(0, 0, 0);
    R.resize(9, 1); 
    R = robotStates[0].block<9, 1>(3, 0); 
    R.resize(3, 3);
    p_world = robotStates[0].block<3, 1>(0, 0) + R * robotStates[1].block<3, 1>(0, 0);

    for (int i = 0 ; i < 1000; ++i) 
    
    {

        // timer.tic();
        simulateStep(tau1);
        p[0] += J_q[0].block(0, 0, 3, num_tendons)*tau1*dt; 
        p[1] += J_q[1].block(0, 0, 3, num_tendons)*tau1*dt; 
        p_world += J_world[1].block(0, 0, 3, num_tendons)*tau1*dt;
        // timer.toc();
        

    }

    robotStates[0] = integrateStates(0);
    robotStates[1] = integrateStates(1);
    std::cout << "robotStates at 0" << ": " << robotStates[0] << "\n\n\n\n";
    std::cout << "robotStates at 1" << ": " << robotStates[1] << "\n\n\n\n";
    std::cout << "p_world:  " << p_world << "\n\n";
    R.resize(9, 1); 
    R = robotStates[0].block<9, 1>(3, 0); 
    R.resize(3, 3);
    p_world = robotStates[0].block<3, 1>(0, 0) + R * robotStates[1].block<3, 1>(0, 0);
    std::cout << p[0] << "\n\n"; 
    std::cout << p[1] << "\n\n"; 
    std::cout << "p_world:  " << p_world << "\n\n";
    

    std::cout << getAffine3d_Distal(0).matrix() << "\n\n";
    std::cout << MathUtils::adjointTransformation(getAffine3d_Distal(0)) << "\n\n\n";

}

Eigen::Affine3d MultistageTDCR_Solver::getAffine3d_Distal(unsigned int stage_num) 

{

    T_Affine3d.translation() = robotStates[stage_num].block<3, 1>(0, 0);
    R.resize(9, 1); 
    R = robotStates[stage_num].block<9, 1>(3, 0);
    R.resize(3, 3);
    T_Affine3d.linear() = R;

    return T_Affine3d; 

}

Eigen::Affine3d MultistageTDCR_Solver::getAffine3d_Proximal(unsigned int stage_num) 

{

    T_Affine3d.translation() = initialConditions[stage_num].block<3, 1>(0, 0);
    R.resize(9, 1); 
    R = initialConditions[stage_num].block<9, 1>(3, 0);
    R.resize(3, 3);
    T_Affine3d.linear() = R;

    return T_Affine3d; 

}

void MultistageTDCR_Solver::setInitialConditions(unsigned int stage_num, Eigen::Matrix<double, 6, 1> ic_force_moment) {

    initialConditions[stage_num].block<6,1>(num_p + num_R, 0) = ic_force_moment;

} 

MultistageTDCR_Solver::~MultistageTDCR_Solver() {};

