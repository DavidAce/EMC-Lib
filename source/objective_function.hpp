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
#include <type_traits>


using namespace Eigen;
class objective_function{
private:
public:
    template <typename func_type,typename boundaryType, typename initial_type,typename... aux_args>
    objective_function(func_type    function_to_minimize,
                       boundaryType lower_boundary,
                       boundaryType upper_boundary,
                       double minimum_tolerance,
                       initial_type initial_conditions,
                       aux_args... optional_data)
            : tolerance(minimum_tolerance),
              aux({optional_data...})
    {
        provided_function = function_to_minimize;
        num_parameters = (int) lower_boundary.size();
        if (num_parameters!= upper_boundary.size()){
            std::cout << "Boundary size mismatch! Exiting" << endl;
            exit(1);
        }
        lower_bound.resize(lower_boundary.size());
        upper_bound.resize(upper_boundary.size());

        for (int i = 0; i < num_parameters; i++){
            lower_bound(i) = lower_boundary(i);
            upper_bound(i) = upper_boundary(i);
            if (lower_bound (i) == upper_bound(i)){
                cout << "Upper and lower boundary [" << i << "] are equal. Exiting" << endl;
                exit(1);
            }
        }

        if (initial_conditions.size() != num_parameters){
            initial_conditions_passed = false;
        }else{
            initial_conditions_passed = true;
            initial.resize(initial_conditions.size());
            for(int i = 0; i < initial_conditions.size();i++){
                initial(i) = initial_conditions(i);
            }
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
    bool initial_conditions_passed;
    ArrayXd initial;

    std::vector<Array<double,Dynamic,Dynamic>>  aux; //Auxiliary data or vectors for use in the fitness function

};


#endif //OBJECTIVE_FUNCTION_H
