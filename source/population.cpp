
#include "population.hpp"


using namespace std;
using namespace EMC_constants;
using namespace Eigen;


std::ostream &operator<<(std::ostream &os, population const &pop) {
	for (int i = 0; i < N; i++) {
		os << pop.guys[i].H << endl;
	}
	return os;
}

void population::getFitness(personality& guy){
	guy.genome.update_parameters();
	guy.H = obj_fun.fitness(guy.genome.parameters);
}


void population::getFitness(const ArrayXd &p, personality &guy){
    //Set all parameters at once with an ArrayXd
    guy.genome.set_parameters(p);
    guy.H = obj_fun.fitness(guy.genome.parameters);
}



void population::find_lowest_guy(){
    lowest_H = guys[0].H;
    lowest_i = 0;
    for (unsigned int i = 0; i < guys.size(); i++){
        if (guys[i].H < lowest_H){
            lowest_H = guys[i].H;
            lowest_i = i;
        }
    }
}

void population::find_lowest_guy_excluding(std::set<unsigned int> &excluded){
    lowest_H = guys[0].H;
    lowest_i = 0;
    for (unsigned int i = 0; i < guys.size(); i++){
        if (excluded.find(i) != excluded.end()) { continue; }
        if (guys[i].H < lowest_H){
            lowest_H = guys[i].H;
            lowest_i = i;
        }
    }
}

void population::insertguy(int from, int to) {
    //Push out the worst bestguy and move everybody up until "to"
    for (int i = 0; i < to; i++) {
        copy(bestguys[i], bestguys[i + 1]);
    }
    //Insert the exceptional guy into the illustrious group of best guys
    copy(bestguys[to], guys[from]);
}

void population::find_elite() {
    //Here we attempt to find someone better in guys that is better than any guy in best_guys
    std::set<unsigned int> tried_guys;
    int j = 0; //Counts how many have been tried
    while (j < N_best) {
        find_lowest_guy_excluding(tried_guys);
        tried_guys.insert(lowest_i);
        j++;
        //By now we should have a winner, check if he's better than any elite
        for (unsigned int i = N_best ;  i-- > 0;) {
            if (lowest_H < bestguys[i].H) {
                insertguy(lowest_i, i);
                break;
            }
            if (lowest_H == bestguys[i].H) {
                break;
            }
        }

    }
}


void population::roulette_select_two_guys(double &Z) {
    //Selected 0 is the guy with good fitness (lowest), selected 1 is random
    Array<double, N, 1> roulette;
    double lucky_number;
    //Start by finding the best guy in the list
    find_lowest_guy();
    Z = 0;
    for (int i = 0; i < N; i++) {
        Z += exp(-(guys[i].H - lowest_H));
        roulette(i) = Z; //Cumulative sum - Low fitness gives large area on roulette
    }
    lucky_number = uniform_double( (double) 0.0 , Z);
    for (int i = 0; i < N; i++) {
        if (lucky_number <= roulette(i)) {
            selected(0) = i;
            break;
        }
    }
    selected(1) = selected(0);
    while (selected(1) == selected(0)) {
        selected(1) = uniform_integer( 0, N - 1);
    }
}

void population::roulette_select_two_guys() {
    //Selected 0 is the guy with good fitness (lowest), selected 1 is random
    Array<double, N, 1> roulette;
    double lucky_number;
    //Start by finding the best guy in the list
    find_lowest_guy();
    double Z = 0;
    for (int i = 0; i < N; i++) {
        Z += exp(-(guys[i].H - lowest_H));
        roulette(i) = Z; //Cumulative sum - Low fitness gives large area on roulette
    }
    lucky_number = uniform_double( (double) 0.0 , Z);
    for (int i = 0; i < N; i++) {
        if (lucky_number <= roulette(i)) {
            selected(0) = i;
            break;
        }
    }
    selected(1) = selected(0);
    while (selected(1) == selected(0)) {
        selected(1) = uniform_integer( 0, N - 1);
    }
}


int population::roulette_select_one_guy() {
    //Selected 0 is the guy with good fitness (lowest), selected 1 is random
    Array<double, N, 1> roulette;
    double lucky_number;
    //Start by finding the best guy in the list
    find_lowest_guy();
    double Z = 0;
    for (int i = 0; i < N; i++) {
        Z += exp(-(guys[i].H - lowest_H));
        roulette(i) = Z; //Cumulative sum - Low fitness gives large area on roulette
    }
    lucky_number = uniform_double( (double) 0.0 , Z);
    for (int i = 0; i < N; i++) {
        if (lucky_number <= roulette(i)) {
            return i;
        }
    }
}


void population::wakeUpNewGuys(){
    for (int i = 0; i < N; i++) {
        copy(newguys[i], guys[i]);
    }
}

void population::wakeUpSnooker(){
    int j = 0;
    for (int i = 0; i < r_num; i++) {
        copy(snookerGuys[i], guys[j++]);
        if (j == N){j = 0;};
    }
}

void population::wakeUpGuys() {
	ArrayXd T(N);
	//Initialize some temperature ladder, here logarithmic.
	T = LogSpaced(N,Tmax, Tmin);
    for (int i = 0; i < N; i++) {
		guys[i].t = T(i); 		//Assign temperatures produced above
        if (obj_fun.initial_conditions_passed){
            guys[i].genome.copy_initial_conditions();
        }else{
            guys[i].genome.randomize_dna();
        }
        guys[i].H = obj_fun.fitness(guys[i].genome.parameters);
//        cout << setprecision(20) << guys[i].H << endl;
    }
}


void population::wakeUpBest() {
//	ArrayXi copied_guys(N_best);
    std::set<unsigned int> copied_guys;
    int j = 0;

    double  lowest_H;
	unsigned int lowest_i; //The position i on the temperature ladder of the guy with lowest H
    while (copied_guys.size() < N_best && j < N_best) {
        lowest_H = guys[0].H;
        lowest_i = 0;
        //Find the best guy yet among guys
        for (unsigned int i = 0; i < N; i++) {
            //Check if i is in skip-list
            if (copied_guys.find(i) != copied_guys.end()){continue;}
            if (guys[i].H < lowest_H) {
                lowest_H = guys[i].H;
                lowest_i = i;
            }
        }
        //By now we should have a winner, copy him to bestguys and add him to skiplist
        copy(bestguys[N_best - j - 1], guys[lowest_i]);
//        cout << "Bestguy: "<< bestguys[N_best - j - 1].H << endl;
        copied_guys.insert(lowest_i);
        j++;
    }

}

ArrayXd population::local_champion_parameters(){
    double champion_H = bestguys[0].H;
    int    champion_i = 0;
    for (unsigned int i = 0; i < bestguys.size(); i++){
        if (bestguys[i].H < champion_H){
            champion_H = bestguys[i].H;
            champion_i = i;
        }
    }
    return bestguys[champion_i].genome.parameters;

}
double population::local_champion_fitness(){
    double champion_H = bestguys[0].H;
    int    champion_i = 0;
    for (unsigned int i = 0; i < bestguys.size(); i++){
        if (bestguys[i].H < champion_H){
            champion_H = bestguys[i].H;
            champion_i = i;
        }
    }
    return bestguys[champion_i].H;

}

int population::local_champion_index(){
    double champion_H = bestguys[0].H;
    int    champion_i = 0;
    for (unsigned int i = 0; i < bestguys.size(); i++){
        if (bestguys[i].H < champion_H){
            champion_H = bestguys[i].H;
            champion_i = i;
        }
    }
    return champion_i;

}



void population::copy(personality &destination, const personality &source) {
	destination.H					= source.H;
	destination.t					= source.t;
	destination.genome.parameters	= source.genome.parameters;
	destination.genome.chromosomes  = source.genome.chromosomes;

	// std::copy(destination.genome.chromosomes, destination.genome.chromosomes + nGenes, source.genome.chromosomes);
}

void population::copy(DNA &destination, const DNA &source) {
	destination.chromosomes  = source.chromosomes;
}


