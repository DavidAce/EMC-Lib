#ifndef SPECIES_H   // if x.h hasn't been included yet...
#define SPECIES_H   //  #define this so the compiler knows it has been included
#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iomanip>
#include <string>
#include <fstream>
#include <time.h>

#include "personality.hpp"
#include "population.hpp"
#include "datafiles.hpp"
#include "constants.hpp"
#include "mymath.hpp"
#include "DNA.hpp"
#include "evolution.hpp"
#include "objective_function.hpp"

using namespace std::chrono;

typedef std::chrono::high_resolution_clock Clock;
class counters {
public:
	int store_counter;
	int generation;
	double simulation_time;
	high_resolution_clock::time_point simulation_tic;
	high_resolution_clock::time_point simulation_toc;
	double evolution_time;
	high_resolution_clock::time_point evolution_tic;
	high_resolution_clock::time_point evolution_toc;
};



class species{
private:
	void zero_counters();
public:

    species(std::function<ArrayXd(ArrayXXd &dos, ArrayXd &E, ArrayXd &M, ArrayXd &T)> func, const int & parameters, const ArrayXd &min_bound, const ArrayXd &max_bound) ;
	//inData obf;				//Experimental data for minimization
    objective_function obj_fun;
    population pop[M];		//Array of separate populations to evolve independently
	outData out;			//Experimental data for minimization
	counters count;
	double champion_fitness();
	double champion_value();
	int champion_number();
	void print_progress();
	void copy(personality &, personality &);
	void wakeUpAll();
};


#endif
