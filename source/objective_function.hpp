//
// Created by david on 2016-09-05.
//

#ifndef OBJECTIVE_FUNCTION_H
#define OBJECTIVE_FUNCTION_H

#include <memory>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "constants.hpp"
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>

using namespace Eigen;
class objective_function{
private:
public:
    template <typename func_type,typename boundaryType, typename... aux_args>
    objective_function(func_type func, boundaryType lower, boundaryType upper, double tol,aux_args... in)
            : tolerance(tol),
              aux({in...})
    {
        provided_function = func;
        num_parameters = (int) lower.size();
        if (num_parameters!= upper.size()){
            std::cout << "Boundary size mismatch! Exiting" << endl;
            exit(1);
        }
        lower_bound.resize(lower.size());
        upper_bound.resize(upper.size());

        for (int i = 0; i < num_parameters; i++){
            lower_bound(i) = lower(i);
            upper_bound(i) = upper(i);
        }

        id = -1;
        name = "";
    };



    typedef std::function<double(objective_function &, ArrayXd &)> providedType ;
    providedType provided_function;

    double fitness(ArrayXd &parameters){
        return provided_function(*this, parameters);
    }

    int id; // Optional: An id for use in parallel (MPI) computations, when this lib is used elsewhere.
    string name;
    ArrayXd lower_bound; //Minimum allowed values for each fitting parameter
    ArrayXd upper_bound; //Maximum allowed values for each fitting parameter
    double tolerance;                           //If the fitness saturates within tolerance, the program terminates
    ArrayXd optimum;
    int num_parameters;
    std::vector<Array<double,Dynamic,Dynamic>>  aux; //Auxiliary data or vectors for use in the fitness function

};


#endif //OBJECTIVE_FUNCTION_H
