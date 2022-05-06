#include <Eigen/Dense>
#include <Eigen/Core>
#include <vector>
#include <iostream>

int main() 

{

    std::vector<Eigen::VectorXi> test_stages;
    std::vector<Eigen::VectorXi> test_stages_Jacobians;

    Eigen::VectorXi stage; 
    stage.resize(2); 
    stage << 0, 1; 
    test_stages.push_back(stage); 
    stage.resize(2); 
    stage << 2, 3; 
    test_stages.push_back(stage); 
    stage.resize(2); 
    stage << 4, 5;
    test_stages.push_back(stage); 
    int vector_size = 0;
    int index = 0; 

    for (int i = 0; i < 3; ++i) 

    {

        for (int k = 3 - i - 1; k >= 0; --k)

        {

            vector_size += test_stages[k].size();

        }

        stage.resize(vector_size);

        for (int ii = 0 + i; ii < 3; ++ii) 

        {

            stage.segment(index, test_stages[ii].size()) = test_stages[ii];
                  
            index += test_stages[ii].size();

        }

        test_stages_Jacobians.push_back(stage);
        vector_size = 0; 
        index = 0;

    }

    for (int i = 0; i < 3; i++) 
    
    {

        std::cout << "Mapped test_stages_Jacobians: " << test_stages_Jacobians[i] << "\n\n";

    }

    


    return 0; 

} 