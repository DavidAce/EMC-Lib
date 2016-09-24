// PID.cpp : Defines the entry point for the console application.
//
#include "EMC.h"
#include <mpi.h>

//boundaries bounds;
void minimize(objective_function & obj_fun){
    //Start MPI
    int mpi_on;
    MPI_Initialized(&mpi_on);
    MPI_Comm MPI_COMM_EMC;
    if (!mpi_on){
        //Multithreaded
        MPI_Init(NULL, NULL);
        MPI_Comm_dup(MPI_COMM_WORLD,  &MPI_COMM_EMC);
    }else{
        //Single threaded
        MPI_Comm_dup(MPI_COMM_SELF,   &MPI_COMM_EMC);
    }
    //Start some constants;
    EMC_constants::nGenes       = obj_fun.num_parameters;
    EMC_constants::geneLength   = maxbits;
    EMC_constants::genomeLength = EMC_constants::nGenes * EMC_constants::geneLength;

    //Initialize populations
    population pop(obj_fun);
    MPI_Comm_dup(MPI_COMM_EMC,  &pop.MPI_COMM_EMC);
    MPI_Comm_rank(MPI_COMM_EMC, &pop.world_ID);           //Establish thread number of this worker
    MPI_Comm_size(MPI_COMM_EMC, &pop.world_size);         //Get total number of threads
    pop.obj_fun.tolerance         = fmax(1e-18, obj_fun.tolerance);


    //Initialize global hall of fame class that contains best champions ever
    hall_of_fame champions;


//    Eigen::initParallel();
	//Start algorithm
	pop.count.simulation_tic = high_resolution_clock::now();
    rng.seed((unsigned long int) pop.world_ID);
    pop.count.timer_store_best += EMC_constants::rate_store_best + 1;
    pop.count.timer_print_progress += EMC_constants::rate_print +1;
    champions.store_champion(pop);
    champions.print_progress(pop);
    while (pop.count.generation < EMC_constants::max_generations &&  !champions.converged(pop)) {
//        cout << "ID: " << pop.world_ID << " Migration" << " Generation = " << pop.count.generation<< endl;

        migration(pop);

//        cout << "ID: " << pop.world_ID << " evolution"<< " Generation = " << pop.count.generation << endl;
        evolve(pop); 			//Evolve all the populations
//        cout << "ID: " << pop.world_ID << " Storing"<<  " Generation = " << pop.count.generation << endl;
        champions.store_champion(pop);
//        cout << "ID: " << pop.world_ID << " Printing"<< " Generation = " << pop.count.generation << endl;
        champions.print_progress(pop);
//        cout << pop.guys[N-1].H << endl;
	}
    //Print final results
    pop.count.timer_print_progress += EMC_constants::rate_print +1;
    champions.print_progress(pop);
    champions.print_final_result(pop);
    //Store final result in objective function
    obj_fun.optimum = champions.champion_parameters;
	//Print timing to console
	pop.count.simulation_toc = high_resolution_clock::now();
	printf("Total time:		%.3f seconds\n", std::chrono::duration<double>(pop.count.simulation_toc - pop.count.simulation_tic).count());
    if(!mpi_on){
        MPI_Finalize();
    }

}