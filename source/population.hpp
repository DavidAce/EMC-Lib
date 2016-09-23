
#ifndef POPULATION_H   // if x.h hasn't been included yet...
#define POPULATION_H   //  #define this so the compiler knows it has been included
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <set>
#include <chrono>
#include <iomanip>
#include <string>
#include <fstream>
#include <time.h>

#include "personality.hpp"
#include "objective_function.hpp"
#include "constants.hpp"
#include "mymath.hpp"
using namespace std;
using namespace EMC_constants;
using namespace Eigen;
//class inData;		//Forward declaration

class paramLine{
	private:
        objective_function &obj_fun;
        ArrayXd  o; //
        ArrayXd  d;
        ArrayXd distance_to_boundaries;

public:
		paramLine(objective_function &ref)
                :obj_fun(ref){
            o.resize(nGenes);
            d.resize(nGenes);
            distance_to_boundaries.resize(2*nGenes);

        };
        double line_min;
        double line_max;
		//Always use "Through" first, to set "o" and "d".
		void Through(ArrayXd &p1, ArrayXd  &p2){ //Points in parameter space p1 and p2
            o = p1;
            d = p2-p1;
            distance_to_boundaries << (obj_fun.upper_bound - o).cwiseQuotient(d), (obj_fun.lower_bound - o).cwiseQuotient(d);
//			distance_to_boundaries = (distance_to_boundaries == distance_to_boundaries).select(distance_to_boundaries, std::numeric_limits<double>::infinity());

			line_max = (distance_to_boundaries < 0).select(distance_to_boundaries.maxCoeff() + 1,  distance_to_boundaries).minCoeff();
            line_min = (distance_to_boundaries > 0).select(distance_to_boundaries.minCoeff() - 1,  distance_to_boundaries).maxCoeff();
//            if( std::isnan(line_max)){
//                cout << p1.transpose() << endl;
//                cout << p2.transpose() << endl;
//                cout << d.transpose() << endl;
//                cout << "line_min = "    << line_min << " line_max = " << line_max << endl;
//                cout << "r is nan!! " << endl;
//                exit(1);
//            }
//            if( line_min * line_max > 0){
//                cout << p1.transpose() << endl;
//                cout << p2.transpose() << endl;
//                cout << d.transpose() << endl;
//                cout << "line_min = "    << line_min << " line_max = " << line_max << endl;
//                cout << "r have the same sign!!! " << endl;
//                exit(1);
//                return;
//            }
//            if (line_min == line_max){
//                cout << p1.transpose() << endl;
//                cout << p2.transpose() << endl;
//                cout << d.transpose() << endl;
//                cout << "line_min = "    << line_min << " line_max = " << line_max << endl;
//                cout << "r are equal!!! " << endl << endl;
//                exit(1);
//                return;
//            }
//
//            if (fmin(line_min,line_max) > 0.0 ){
//                cout << p1.transpose() << endl;
//                cout << p2.transpose() << endl;
//                cout << d.transpose() << endl;
//                cout << "line_min = "    << line_min << " line_max = " << line_max << endl;
//                cout << "line_min too large!!! " << endl << endl;
//                exit(1);
//            }
//
//            if (fmax(line_min,line_max) < 1.0 ){
//                cout << p1.transpose() << endl;
//                cout << p2.transpose() << endl;
//                cout << d.transpose() << endl;
//
//                cout << "line_min = "    << line_min << " line_max = " << line_max << endl;
//                cout << "line_max too too small!!! "<< endl << endl;
//                exit(1);
//
//
//            }
//            if (line_min == 0 || line_max == 0){
//                cout << p1.transpose() << endl;
//                cout << p2.transpose() << endl;
//                cout << d.transpose() << endl;
//                cout << "line_min = "    << line_min << " line_max = " << line_max << endl;
//
//                cout << "r are zero!!! " << endl;
//
//                exit(1);
//
//            }

        }
        ArrayXd  pointAt(const double r){
            return o + d*r;
		}
//		double line_max() {
//            return (obj_fun.upper_bound - o).cwiseQuotient(d).minCoeff();
//        }
//	    double line_min() { //Get the largest r within upper boundary
//		    return -(obj_fun.lower_bound - o).cwiseQuotient(-d).minCoeff();
//	    }
};


typedef std::chrono::high_resolution_clock Clock;
class counters {
public:
    counters(){
        generation          = 0;
        timer_migration     = 0;
        timer_store_best    = 0;
        timer_print_progress= 0;
        timer_check_conv    = 0;
        simulation_time     = 0;
        evolution_time      = 0;
    }
	int store_counter;
	int store_last_since;
	int generation;
    int timer_migration;
    int timer_store_best;
    int timer_check_conv;
    int timer_print_progress;
	double simulation_time;
	Clock::time_point simulation_tic;
	Clock::time_point simulation_toc;
	double evolution_time;
	Clock::time_point evolution_tic;
	Clock::time_point evolution_toc;
};




class population{
private:

	void wakeUpGuys();
	void wakeUpNewGuys();
	void wakeUpBest();
public:
	population(objective_function &ref)
			:obj_fun        (ref),
			 guys 	        (N,ref),
			 newguys        (N,ref),
			 bestguys       (N_best,ref),
			 snookerGuys    (r_num,ref),
             line           (ref)
    {
                wakeUpGuys();
                wakeUpNewGuys();
                wakeUpBest();
    };
    objective_function &obj_fun;
    vector<personality> guys; 			     //Make an array of N guys
	vector<personality> newguys; 			 //Make a temporary array of N guinneapigs
	vector<personality> bestguys;            //Make an array of N/10 good performers
	vector<personality> snookerGuys;         // (r_num, personality(true));//Make an array of r_num snooker guys

    paramLine line;	//for doing the snooker crossover
	double lowest_H;
    unsigned int    lowest_i;
    Array2i selected;

    counters count;
    int world_ID;     //Numbers for MPI communication
    int world_size;   //Numbers for MPI communication
    void getFitness(personality &guy);
    void getFitness(const ArrayXd  &point, personality &guy);
    void copy(personality &destination, const personality & source);
    void copy(DNA& destination, const DNA& source);
    void find_lowest_guy();

    void find_lowest_guy_excluding(std::set<unsigned int> &excluded);
    void insertguy(int from, int to);
    void find_elite();
    void roulette_select_two_guys(double &Z);
    void roulette_select_two_guys();
    int  roulette_select_one_guy();
    ArrayXd local_champion_parameters();
    double  local_champion_fitness();
    int     local_champion_index();
	int operator()() { //Return the bit at a.
		return 0;
	}
	friend ostream &operator<<(std::ostream &os, population const &);


};

#endif
