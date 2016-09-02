//
// Created by david on 9/2/16.
//

#include "objective_function.hpp"


objective_function::objective_function(std::function<ArrayXd(ArrayXXd &dos, ArrayXd &E, ArrayXd &M, ArrayXd &T)> func_,
                                       const int & parameters_,
                                       const ArrayXd &min_bound_,
                                       const ArrayXd &max_bound_)
:the_function(func_),
 parameters(parameters_),
 min_bound(min_bound_),
 max_bound(max_bound_)
{


}
