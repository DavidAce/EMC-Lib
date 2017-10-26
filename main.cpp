//
// Created by david on 9/2/16.
//

#include "EMC.h"
#include <mpi.h>
#include <unsupported/Eigen/CXX11/Tensor>

double my_example_function(objective_function &obj_fun, ArrayXd &inputParameters){
    //ArrayXd has the same size as your lower and upper boundaries. You can map it to whatever using the following:
    //     auto my_tensor = TensorMap<Tensor<double,3>>(inputParameters.data(), 3,3,3);  //(3 is rank, and 3,3,3 the sizes of the tensor)
    //     auto my_Matrix = Map<MatrixXd>(inputParameters.data(), 3,3);  // 3,3 is the rows and columns of the matrix

    double accumulator = 0;
    for (int i = 0; i < inputParameters.size() ; i++){
        accumulator += inputParameters(i)*inputParameters(i);
    }
    //Notice below how obj_fun.aux contains the exampleData you passed earlier.
    return sqrt(accumulator) + obj_fun.aux[0].sum() - obj_fun.aux[1].sum(); //Minimum at 0
};


int main(){
    cout << "This example minimizes a 2D paraboloid f(x,y) = sqrt(x^2 + y^2) with solution (0,0)" << endl;

    //This example data is there to showcase the syntax of the program.
    Eigen::ArrayXd exampleData0(4);
    Eigen::ArrayXd exampleData1(4);

    exampleData0 << 1,2,4,3;
    exampleData1 << 3,2,4,1;

    //Mandatory arrays! Types Eigen types (Tensor, Array,Matrix,) or std::vector<double>, etc.
    //just make sure you can do .size() on it.
    Eigen::Tensor<double,3> lower_bound(3,3,3);
    Eigen::Tensor<double,3> upper_bound(3,3,3);
    Eigen::Tensor<double,3> initial_condition(3,3,3);

    lower_bound.setConstant(-10);
    upper_bound.setConstant(10) ;
    initial_condition.setConstant(1) ;
    double tolerance = 1e-12;   //The program terminates once the fitness does not improve beyond this tolerance


    //This function takes arguments lower_bound, upper_bound, tolerance,
    //and then any number of auxiliary double arrays (ArrayXd or ArrayXXd) of any size.
    //These can be accessed from obj_fun.aux[0], obj_fun.aux[1] ... etc, in the same order as given
    objective_function obj_fun(my_example_function,lower_bound, upper_bound,tolerance , initial_condition,exampleData0,exampleData1);
    obj_fun.id = 0;
    obj_fun.name = "Parabola | ";
    //Next we pass a function to minimize.
    // Its output is a double (the fitness)
    // Its input  is a reference objective_function &obj_fun, and an ArrayXd &parameters with  fitting parameters
    minimize(obj_fun);
    return 0;
}
