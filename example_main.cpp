//
// Created by david on 9/2/16.
//

#include <Eigen/Dense>
#include <Eigen/Core>
#include "source/EMC.h"
#include "objective_function.hpp"

using namespace Eigen;

long double my_example_function(objective_function &obj_fun, Array<long double, Dynamic, 1> &inputParameters){
    return sqrt(inputParameters.cwiseAbs2().sum()) + obj_fun.aux[0].sum() - obj_fun.aux[1].sum(); //Minimum at 0
};



int main(){
    cout << "This example minimizes a 2D paraboloid f(x,y) = sqrt(x^2 + y^2) with solution (0,0)" << endl;

    //This example data is there to showcase the syntax of the program.
    ArrayXd exampleData0(3);
    ArrayXd exampleData1(3);

    exampleData0 << 1,2,3;
    exampleData1 << 3,2,1;

    //Mandatory arrays!
    ArrayXd lower_bound(2);
    ArrayXd upper_bound(2);

    lower_bound << -10, -10;
    upper_bound << 10, 10;
    //The program terminates once the fitness does not improve beyond this tolerance
    double tolerance = 1e-16;

    //This function takes arguments lower_bound, upper_bound, tolerance,
    //and then any number of auxiliary double arrays (ArrayXd or ArrayXXd) of any size.
    //These can be accessed from obj_fun.aux[0], obj_fun.aux[1] ... etc, in the same order as given
    objective_function obj_fun(lower_bound, upper_bound,tolerance ,exampleData0,exampleData1);

    //Next we pass a function to minimize.
    // Its output is a double (the fitness)
    // Its input  is a reference objective_function &obj_fun, and an ArrayXd &parameters with  fitting parameters
    obj_fun.provided_function = my_example_function;
    minimize(obj_fun);
    return 0;
}