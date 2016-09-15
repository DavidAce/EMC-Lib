// PID.cpp : Defines the entry point for the console application.
//
#include "EMC.h"


//boundaries bounds;
void minimize(objective_function & obj_fun){
	//Start up some files and folder for saving out
    obj_fun.tolerance           = fmax(1e-18, obj_fun.tolerance);
    EMC_constants::nGenes 	   	= obj_fun.parameters;
    EMC_constants::geneLength   = 2+min(58,(int)ceil(-log(obj_fun.tolerance)/log(2)));
    EMC_constants::genomeLength = EMC_constants::nGenes * EMC_constants::geneLength;
	if (obj_fun.threads >= 0){
        EMC_constants::M = obj_fun.threads;
        omp_set_num_threads(obj_fun.threads);
    }else{
        EMC_constants::M = omp_get_num_threads();
    }
    cout << "OpenMP Threads: " << EMC_constants::M<< endl;

    Eigen::initParallel();
    species sp(obj_fun);

    rng.seed(EMC_constants::seed);
	//Start algorithm
	sp.count.simulation_tic = high_resolution_clock::now();
	// sp.count.evolution_tic = clock();				//Start timer
	#pragma omp parallel
    rng.seed(EMC_constants::seed + (unsigned long)omp_get_thread_num());
	while (sp.count.generation < EMC_constants::max_generations &&  !sp.below_tolerance()) {
		#pragma omp single nowait
		{

            sp.print_progress();
            sp.store_best_fitness();
        if (EMC_constants::M > 1) {
            if (uniform_double_1() < qmig) {
                migration(sp);
            }
        }
		}
		#pragma omp for nowait
		for (int i = 0; i < M; i++) {
			evolve(sp.pop[i]); 			//Evolve all the populations
		}
		

		
	}
	//Print final parameters
    sp.print_progress(true);
    if (obj_fun.id >= 0){cout << "ID: " << obj_fun.id << " ";}
    cout << "Best Parameters: "
		 << sp.pop[sp.champion_number()].bestguys[N_best - 1].genome.parameters.transpose() <<"     | " ;
	//Print timing to console
	sp.count.simulation_toc = high_resolution_clock::now();
	printf("Total time:		%.3f seconds\n", std::chrono::duration<double>(sp.count.simulation_toc - sp.count.simulation_tic).count());
	obj_fun.optimum = sp.pop[sp.champion_number()].bestguys[N_best - 1].genome.parameters ;
}