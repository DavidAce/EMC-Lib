//
// Created by david on 9/2/16.
//

#ifndef EMC_CLION_OBJECTIVE_FUNCTION_H
#define EMC_CLION_OBJECTIVE_FUNCTION_H
#include <Eigen/Dense>
#include <Eigen/Core>
#include "personality.hpp"
using namespace std;
using namespace Eigen;

class objective_function {
private:
    std::function<ArrayXd(ArrayXXd &dos, ArrayXd &E, ArrayXd &M, ArrayXd &T)> the_function;

public:

    objective_function(std::function<ArrayXd(ArrayXXd &dos, ArrayXd &E, ArrayXd &M, ArrayXd &T)> func, const int & parameters, const ArrayXd &min_bound, const ArrayXd &max_bound);

    const int parameters;
    const ArrayXd min_bound;
    const ArrayXd max_bound;

    double fitnessTest(personality &guy) {
//        cout << "Hej 2 1 1" << endl;

        double H = 0
        ;
        guy.value = H; //Record the total distance between mappings
        //(DO NOT CHANGE) I propose to use the function below to make H small overall, and sharp close to H = 0


        H = -1 / log(H + log_param) + log_const + pow(log( 1/(H + 1)), 2);
        return H;
    }

};


#endif //EMC_CLION_OBJECTIVE_FUNCTION_H
