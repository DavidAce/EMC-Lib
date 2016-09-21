//
// Created by david on 9/2/16.
//

#include "EMC.h"
#include <mpi.h>
#include <unsupported/Eigen/CXX11/Tensor>

double my_example_function(objective_function &obj_fun, Eigen::Tensor<double, 3> &inputParameters){
    double accumulator = 0;
    for (int i = 0; i < inputParameters.size() ; i++){
        accumulator += inputParameters(i)*inputParameters(i);
    }
    return sqrt(accumulator) + obj_fun.aux[0].sum() - obj_fun.aux[1].sum(); //Minimum at 0
};


int main(){
    cout << "This example minimizes a 2D paraboloid f(x,y) = sqrt(x^2 + y^2) with solution (0,0)" << endl;
    MPI_Init(NULL, NULL);
    int world_ID,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);           //Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);         //Get total number of threads


    //This example data is there to showcase the syntax of the program.
    Eigen::ArrayXd exampleData0(4);
    Eigen::ArrayXd exampleData1(4);

    exampleData0 << 1,2,4,3;
    exampleData1 << 3,2,4,1;

    //Mandatory arrays!
    Eigen::Tensor<double,3> lower_bound(3,3,3);
    Eigen::Tensor<double,3> upper_bound(3,3,3);

    lower_bound.setConstant(-10);
    upper_bound.setConstant(10) ;
    //The program terminates once the fitness does not improve beyond this tolerance
    double tolerance = 1e-16;

    //This function takes arguments lower_bound, upper_bound, tolerance,
    //and then any number of auxiliary double arrays (ArrayXd or ArrayXXd) of any size.
    //These can be accessed from obj_fun.aux[0], obj_fun.aux[1] ... etc, in the same order as given
    objective_function obj_fun(my_example_function,lower_bound, upper_bound,tolerance ,exampleData0,exampleData1);
    obj_fun.id = 0;
    obj_fun.name = "Parabola | ";
    obj_fun.threads = 4;
    //Next we pass a function to minimize.
    // Its output is a double (the fitness)
    // Its input  is a reference objective_function &obj_fun, and an ArrayXd &parameters with  fitting parameters
    minimize(obj_fun);
    MPI_Finalize();

    return 0;
}