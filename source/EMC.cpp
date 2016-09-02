// PID.cpp : Defines the entry point for the console application.
//
#include "EMC.h"

RNGType rng;
//boundaries bounds;
ArrayXd minimization(std::function<ArrayXd(ArrayXXd &dos, ArrayXd &E, ArrayXd &M, ArrayXd &T)> func, const int & parameters, const ArrayXd &min_bound, const ArrayXd &max_bound){
	//Start up some files and folder for saving out

    constants::upper_bound = max_bound;
    constants::lower_bound = min_bound;
    constants::nGenes = parameters;
    constants::genomeLength = parameters*constants::geneLength;
    Eigen::initParallel();
    species sp(func, parameters, min_bound, max_bound);
    cout << "Hej 1" << endl;
    cout << sp.pop[0] << endl;
    cout << sp.pop[0].guys[0].genome << endl;

    rng.seed(6);
	//Start algorithm
	sp.count.simulation_tic = high_resolution_clock::now();
	// sp.count.evolution_tic = clock();				//Start timer
	#pragma omp parallel
	while (sp.count.generation < generations &&  sp.champion_fitness() > lowest_H) {
		#pragma omp single nowait
		{

            sp.print_progress();
		if (uniform_double(&rng, 0, 1) < qmig) {
			migration(sp);
		}
		}
		#pragma omp for nowait
		for (int i = 0; i < M; i++) {
			evolve(sp.pop[i], sp.obj_fun); 			//Evolve all the populations
		}
		

		
	}

	//sp.count.evolution_time += (double)(sp.count.evolution_toc - sp.count.evolution_tic) / CLOCKS_PER_SEC; //Collect computation time

	sp.out.print_to_file(sp);
	//Print final parameters
	cout << endl << "Best Parameters: " 
		 << sp.pop[sp.champion_number()].bestguys[N_best - 1].genome.parameters.transpose() << endl;
	//Print timing to console
	sp.count.simulation_toc = high_resolution_clock::now();
	printf("\nTotal time:		%.3f seconds\n", std::chrono::duration<double>(sp.count.simulation_toc - sp.count.simulation_tic).count());
	return sp.pop[sp.champion_number()].bestguys[N_best - 1].genome.parameters.transpose() ;
}