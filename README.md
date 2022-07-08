# **Tendon-driven Continuum Robot based on Cosserat Rod Theory**

## **Description**

## **Installation**

## **Terminology**

## **Minimal Example**

#### Some Python and C++ skills will be required in order to create a custom continuum robot model.

#### 1. Define continuum robot's parameters.

#### The minimum parameters required are: 

    a. Backbone material properties such as Young's Modulus, Shear Modulus.
    b. Length of backbone, and the length of each section/stage.
    c. Mass distribution of the backbone. 
    d. Number of tendons/cables and the positions of each routing in 3D cartesian space. 
    e. Geometry of the backbone, i.e. hollow rod, solid rod. If no geometry is found, then the bending and shear matrix must be provided. 

#### 2. Generated c code using Python script, found under src. 

#### The steps can be broken down into: 

    a. Create a "robot_dict" that requires:
        - Type of rod used, 
        - Geometry, 
        - Material properties.
    b. Create another dictionary that requires :
        - Number of tendons,
        - Routings in 3D space,
        - Stage lengths, 
        - Stage tendons,
          - Stage tendons simply means which tendons affect which stage/section of the continuum robot.
        - Number of integration steps (Generally longer rods need more integration steps for reconstruct its 3D representation.)
    c. Then, we can use the "Robot_Parameter_Builder" and "Robot_Type_Builder" templates to create the c generated code required for the C++ TDCR_Interface.


#### 3. Create the MPC controller object. 

#### 4. As of this version, we need to create `TDCR_Interface`, which requires both `MultistageTDCR_Solver` and `ControllerInterface`.

    a. We first start by including the required header files for the integrators.
    b. Next, we create IntegrationInterface objects for the MultistageTDCR_Solver. 
    c. Next, we define the stage tendons and routings. Note: They should be the same as described in step 2b!
    d. We can now create the MultistageTDCR_Solver, ControllerInterface and TDCR_Interface.
    e. Then, we use solveForwardKinematics and simulateStep to ensure that the robotStates in the MultistageTDCR_Solver object is valid. TODO: Remove the need for this step.
    f. Create the controlStatesMatrix (CSM) and controlled states (CS) matrices/vectors to set the controller's dimensions. 
    e. Run the simulation/control loop.