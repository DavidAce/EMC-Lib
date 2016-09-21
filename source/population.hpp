
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
        Array<double,Dynamic,1>  o;
        Array<double,Dynamic,1>  d;
	public:
		paramLine(objective_function &ref)
                :obj_fun(ref){
            o.resize(nGenes);
            d.resize(nGenes);
        };
		//Always use "Through" first, to set "o" and "d".
		void Through(Tensor<double,3> &p1, Tensor<double,3>  &p2){ //Points in parameter space p1 and p2
            for(int i = 0;i < nGenes; i++){
				o(i) = p1(i);
				d(i) = p2(i) - p1(i);
			}
		}
        template<typename type>
        Array<type,Dynamic,1>  pointAt(const type r){
            Array<type,Dynamic,1>  v(nGenes);
			for(int i = 0;i < nGenes; i++){
				v(i) = o(i)+d(i)*r;
			}
			return v;
		}
        ArrayXd subtract(double t[],  double a[]){
            ArrayXd result(nGenes);
            for (int i = 0; i < nGenes ; i++){
                result (i) = (int) ( t[i] - a[i]);
            }
            return result;
        }
//		double distance(Tensor<double,Dynamic> &p1, Tensor<double,Dynamic>  &p2){
//			return sqrt((p2.cast<double>()-p1.cast<double>()).square().sum());
//		}
		double line_max() { //Get the largest r within upper boundary Bu and lower boundary Bl
            return  subtract(obj_fun.upper_bound.data(), o.data()).cast<double>().cwiseQuotient(d.cast<double>()).minCoeff();
//            return (obj_fun.upper_bound.cast<double>() - o).cwiseQuotient(d).minCoeff();
		}
	    double line_min() { //Get the largest r within upper boundary Bu and lower boundary Bl
		    return  -subtract(obj_fun.lower_bound.data(), o.data()).cast<double>().cwiseQuotient((-d.cast<double>()).eval()).minCoeff();
	    }
};


typedef std::chrono::high_resolution_clock Clock;
class counters {
public:
	int store_counter;
	int store_last_since;
	int generation;
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
	objective_function &obj_fun;
public:
	population(objective_function &ref)
			:obj_fun        (ref),
			 guys 	        (N,ref),
			 newguys        (N,ref),
			 bestguys       (N_best,ref),
			 snookerGuys    (r_num,ref),
             line           (ref),
			 generation(0)
	{
                wakeUpGuys();
                wakeUpNewGuys();
                wakeUpBest();
//		wakeUpPop();
    };
    vector<personality> guys; 			     //Make an array of N guys
	vector<personality> newguys; 			 //Make a temporary array of N guinneapigs
	vector<personality> bestguys;            //Make an array of N/10 good performers
	vector<personality> snookerGuys;         // (r_num, personality(true));//Make an array of r_num snooker guys

    paramLine line;	//for doing the snooker crossover
    int generation;  //Number of generations for this population
	counters count;

    void getFitness(personality &guy);
    void getFitness(const Array<double,Dynamic,1>  &point, personality &guy);
    void getFitness(const Tensor<double,3> &point, personality &guy);

    void copy(personality &destination, const personality & source);
    void copy(DNA& destination, const DNA& source);


	int operator()() { //Return the bit at a.
		return 0;
	}
	friend ostream &operator<<(std::ostream &os, population const &);


};

#endif
