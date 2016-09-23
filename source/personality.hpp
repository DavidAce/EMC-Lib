
#ifndef PERSONALITY_H   // if x.h hasn't been included yet...
#define PERSONALITY_H   //  #define this so the compiler knows it has been included
#include <Eigen/Dense>
#include <Eigen/Core>
#include <bitset>
#include <iostream>
#include <random>
#include <cfloat>
#include "constants.hpp"
#include "DNA.hpp"
#include "objective_function.hpp"

using namespace std;
using namespace EMC_constants;
using namespace Eigen;

class personality {
private:
    objective_function &obj_fun;
public:

	personality(objective_function &ref)
            :obj_fun(ref),
             born(0),
             genome(ref)
	{
        EMC_constants::nGenes = obj_fun.num_parameters;
        EMC_constants::nGenes = obj_fun.num_parameters;
        EMC_constants::geneLength   = 2+min(54,(int)ceil(-log(obj_fun.tolerance)/log(2)));
        EMC_constants::genomeLength = EMC_constants::nGenes * EMC_constants::geneLength;
        genome.chromosomes.resize((unsigned int)nGenes);
        genome.parameters.resize ((unsigned int)nGenes);
        if (obj_fun.initial_conditions_passed){
            genome.copy_initial_conditions();
        }else{
            genome.randomize_dna();

        }
//        H = ref.fitness(genome.parameters); //Set fitness
//        cout << H << endl;

	};
	personality(objective_function &ref, bool toggle) :
            obj_fun(ref),
            born(0),
            genome(ref,toggle) {
	};
	double t; 							//temperature
	int born;							//Generation when DNA first emerged
	DNA genome;							//Contains binary and real representation of parameters
	double H;					//Fitness, or energy
	friend ostream &operator<<(std::ostream &os, personality const &);

};
#endif
