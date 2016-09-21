//
// Created by david on 9/20/16.
//

#ifndef EMC_LIB_PARALLEL_ALGORITHMS_H
#define EMC_LIB_PARALLEL_ALGORITHMS_H

#include "population.hpp"
#include "constants.hpp"
#include "iostream"
#include <mpi.h>

void migration(population &pop);


class hall_of_fame{
private:
    void  global_champion_parameters(population &pop);
    void  global_champion_fitness_and_index(population &pop);
public:
    hall_of_fame(){
        latest_history_diff = 1;
    }
    double          champion_fitness;
    int             champion_pop_index;
    ArrayXd         champion_parameters;
    vector<double>  previous_champions_fitness;
    double          latest_history_diff;



    void update_global_champion(population &pop);
    void store_champion(population &pop);
    void print_progress(population &pop);
    void print_final_result(population &pop);
    bool converged(population &pop);


};



#endif //EMC_LIB_PARALLEL_ALGORITHMS_H
