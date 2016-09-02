//
// Created by david on 9/2/16.
//

#ifndef EMC_CLION_EMC_H
#define EMC_CLION_EMC_H
#include <stdio.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <bitset>
#include <chrono>
#include "constants.hpp"
#include "randomFunctions.hpp"
#include "DNA.hpp"
#include "personality.hpp"
#include "population.hpp"
#include "species.hpp"
#include "datafiles.hpp"
#include "evolution.hpp"
#include "objective_function.hpp"
using namespace Eigen;
using namespace std;
using namespace constants;
using namespace std::chrono;

//Minimization should be a function that takes a pointer to a function,
//An int with the number of parameters to fit, and two arrays with boundaries min and max for each parameter.
ArrayXd minimization(std::function<ArrayXd(ArrayXXd &dos, ArrayXd &E, ArrayXd &M, ArrayXd &T)> func, const int & parameters, const ArrayXd &min_bound, const ArrayXd &max_bound);

#endif //EMC_CLION_EMC_H
