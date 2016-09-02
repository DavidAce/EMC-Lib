
#ifndef EVOLUTION_H
#define EVOLUTION_H  

#include "personality.hpp"
#include "objective_function.hpp"
#include "population.hpp"
#include "species.hpp"

class population;	//Forward declaration
class personality;	//Forward declaration
extern void roulette_select(personality [], int [], double *, double );
extern void bitselector_smartCopy(population &, int [], int []);
extern void mutation(population &, objective_function&);
extern void mutation_elite(population &, objective_function&);
extern void crossover(population &, objective_function&);
extern void crossover_elite(population &, objective_function&);
extern void crossover_smartCopy(population &, objective_function&);
extern void crossover_snooker(population &, objective_function&);
extern void exchange(population &);
extern void migration(species &);
extern void insertguy(population &, int , int );
extern void find_elite(population &);
extern void evolve(population &, objective_function &);

#endif