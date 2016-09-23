
#ifndef DNA_H   // if x.h hasn't been included yet...
#define DNA_H   //  #define this so the compiler knows it has been included
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <bitset>
#include <vector>
#include <iostream>
#include <iomanip>
#include <random>
#include <memory>
#include <climits>
#include <algorithm>
#include "constants.hpp"
#include "randomFunctions.hpp"
#include "objective_function.hpp"

using namespace std;
using namespace EMC_constants;
using namespace EMC_rnd;
using namespace Eigen;
class DNA {
private:
    double bin2dec(const int);
    bitset<maxbits> dec2bin(const int);

//	vector<bool> dec2bin(const int);

//    Array<long double, Dynamic,1>   random_parameters();
//    vector< bitset<maxbits> >       random_chromosomes();



public:
    DNA(objective_function &ref):obj_fun(ref)
    {

    }

    DNA(objective_function &ref, bool ):obj_fun(ref) {
        EMC_constants::nGenes = obj_fun.num_parameters;
        EMC_constants::geneLength   = 2+min(54,(int)ceil(-log(obj_fun.tolerance)/log(2)));
        EMC_constants::genomeLength = EMC_constants::nGenes * EMC_constants::geneLength;
        chromosomes.resize((unsigned int)nGenes);
        parameters.resize ((unsigned int)nGenes);
    }
    objective_function &obj_fun;

    ArrayXd parameters;
//    Array<long double, Dynamic ,1> parameters;						  //Decimal representation
    vector< bitset<maxbits> > chromosomes; //Binary representation

    bool operator== (const DNA& target);
    int operator()(int);
    friend ostream &operator<<(std::ostream &os, DNA const &);
    void flip_loci(const int);
    void flip_loci(ArrayXi &);
    void flip_loci(Ref<ArrayXi> );
    void copy_loci(const int, const int);

    void set_parameter(const int,const double);
    void set_parameters(const  ArrayXd  &p);
    void update_parameters();
    void update_chromosomes();
    void randomize_dna();
    void copy_initial_conditions();
};

#endif
