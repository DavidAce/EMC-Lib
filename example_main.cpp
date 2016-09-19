//
// Created by david on 9/2/16.
//

#include "EMC.h"


long double my_example_function(objective_function &obj_fun, Eigen::Array<long double, Dynamic, 1> &inputParameters){
    return sqrt(inputParameters.cwiseAbs2().sum()) + obj_fun.aux[0].sum() - obj_fun.aux[1].sum() +1e18; //Minimum at 0
};



int main(){
    cout << "This example minimizes a 2D paraboloid f(x,y) = sqrt(x^2 + y^2) with solution (0,0)" << endl;

    //This example data is there to showcase the syntax of the program.
    Eigen::ArrayXd exampleData0(3);
    Eigen::ArrayXd exampleData1(3);

    exampleData0 << 1,2,3;
    exampleData1 << 3,2,1;

    //Mandatory arrays!
    Eigen::ArrayXd lower_bound(2);
    Eigen::ArrayXd upper_bound(2);

    lower_bound << -10, -10;
    upper_bound <<  10, 10;
    //The program terminates once the fitness does not improve beyond this tolerance
    double tolerance = 1e-16;

    //This function takes arguments lower_bound, upper_bound, tolerance,
    //and then any number of auxiliary double arrays (ArrayXd or ArrayXXd) of any size.
    //These can be accessed from obj_fun.aux[0], obj_fun.aux[1] ... etc, in the same order as given
    objective_function obj_fun(my_example_function,lower_bound, upper_bound,tolerance ,exampleData0,exampleData1);
    obj_fun.id = 0;
    obj_fun.threads = 1;
    //Next we pass a function to minimize.
    // Its output is a double (the fitness)
    // Its input  is a reference objective_function &obj_fun, and an ArrayXd &parameters with  fitting parameters
    minimize(obj_fun);
    return 0;
}