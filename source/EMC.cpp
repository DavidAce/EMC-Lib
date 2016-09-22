// PID.cpp : Defines the entry point for the console application.
//
#include "EMC.h"
#include <mpi.h>

//boundaries bounds;
void minimize(objective_function & obj_fun){
    //Start MPI
    MPI_Init(NULL, NULL);
    int world_ID,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_ID);           //Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);         //Get total number of threads
    //Initialize populations
    population pop(obj_fun);
    pop.obj_fun.tolerance           = fmax(1e-18, obj_fun.tolerance);
    pop.world_ID    = world_ID;
    pop.world_size  = world_size;
    //Initialize global hall of fame class that contains best champions ever
    hall_of_fame champions;

    //Start some constants;
    EMC_constants::M            = pop.world_size;
    Eigen::initParallel();
	//Start algorithm
	pop.count.simulation_tic = high_resolution_clock::now();
    rng.seed((unsigned long int) pop.world_ID);
    while (pop.count.generation < EMC_constants::max_generations &&  !champions.converged(pop)) {
//        cout << "ID: " << pop.world_ID << " Migration" << " Generation = " << pop.count.generation<< endl;
        migration(pop);
//        cout << "ID: " << pop.world_ID << " evolution"<< " Generation = " << pop.count.generation << endl;
        evolve(pop); 			//Evolve all the populations
//        cout << "ID: " << pop.world_ID << " Storing"<<  " Generation = " << pop.count.generation << endl;
        champions.store_champion(pop);
//        cout << "ID: " << pop.world_ID << " Printing"<< " Generation = " << pop.count.generation << endl;
        champions.print_progress(pop);
	}
    //Print final results
    champions.print_final_result(pop);
    //Store final result in objective function
    obj_fun.optimum = champions.champion_parameters;
	//Print timing to console
	pop.count.simulation_toc = high_resolution_clock::now();
	printf("Total time:		%.3f seconds\n", std::chrono::duration<double>(pop.count.simulation_toc - pop.count.simulation_tic).count());
    MPI_Finalize();

}