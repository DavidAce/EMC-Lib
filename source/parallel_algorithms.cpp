//
// Created by david on 9/20/16.
//

#include "parallel_algorithms.h"
void migration(population &pop) {
    //First select randomly m out of M pobulations, subject to migration
    //Then select randomly  n (n=1?) out of N guys to migrate
    //Accept the new population based on the new population fitness score
    if (pop.count.timer_migration++ > rate_migration) {
        pop.count.timer_migration = 0;
        int go_ahead = uniform_integer_1(); //Decide with 50% chance whether you migrate this time.
        ArrayXi whos_up_for_migration(M);

        MPI_Allgather(&go_ahead, 1, MPI_INT, whos_up_for_migration.data(), 1, MPI_INT, MPI_COMM_WORLD);
        //Make sure there are at least two who wish to migrate
        if (whos_up_for_migration.sum() < 2) { return; }

        //Now the point is to make a cyclic migration, so it doesnt matter if there are even or odd number
        //of migrators.
        if (go_ahead == 1) {
            int recv_from, send_to;
            int i = pop.world_ID;
            while (true) {
                i++;
                if (i == pop.world_size) { i = 0; }
                if (whos_up_for_migration(i) == 1) {
                    send_to = i;
                    break;
                }
            }
            i = pop.world_ID;
            while (true) {
                i--;
                if (i == -1) { i = pop.world_size - 1; }
                if (whos_up_for_migration(i) == 1) {
                    recv_from = i;
                    break;
                }
            }
            //Now we know exactly who to send/receive to and from
            //Select a guy to send using the roulette wheel.
            int n = pop.roulette_select_one_guy(); //Guy chosen to migrate
            MPI_Sendrecv_replace(pop.newguys[n].genome.parameters.data(), nGenes, MPI_DOUBLE, send_to, pop.world_ID,
                                 recv_from, recv_from, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv_replace(&pop.newguys[n].H, 1, MPI_DOUBLE, send_to, pop.world_ID, recv_from, recv_from,
                                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (uniform_double_1() < exp(-(pop.newguys[n].H - pop.guys[n].H) / pop.newguys[n].t)) {
                //Keep
                pop.newguys[n].born = pop.count.generation;
                pop.newguys[n].genome.update_chromosomes();
                pop.copy(pop.newguys[n], pop.guys[n]);
                pop.guys[n].born = pop.count.generation;


            } else {
                //Revert
                pop.copy(pop.guys[n], pop.newguys[n]);
            }
        }
    }
}



void hall_of_fame::global_champion_parameters(population &pop) {
    champion_parameters = pop.bestguys[N_best - 1].genome.parameters;
    MPI_Bcast(champion_parameters.data(), nGenes, MPI_DOUBLE, champion_pop_index, MPI_COMM_WORLD);
}

//void hall_of_fame::global_champion_fitness(population &pop) {
//    double champion_fitness = pop.bestguys[N_best - 1].H;
//    MPI_Bcast(&champion_fitness, 1, MPI_DOUBLE, champion_pop_index, MPI_COMM_WORLD);
//}

void hall_of_fame::global_champion_fitness_and_index(population &pop) {
//    double local_champion_H = pop.local_champion_fitness();
    struct {
        double champ;
        int    index;
    } in, out;
    in.champ = pop.local_champion_fitness();
    in.index = pop.world_ID;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    champion_fitness = out.champ;
    champion_pop_index  = out.index;
//    MPI_Allreduce(&local_champion_H, &champion_pop_index, 1, MPI_DOUBLE, MPI_MINLOC, MPI_COMM_WORLD);
}

void hall_of_fame::update_global_champion(population &pop) {
    global_champion_fitness_and_index(pop);
    global_champion_parameters(pop);
}


void hall_of_fame::print_progress(population &pop) {
    if (pop.count.timer_print_progress++ > rate_print) {
        pop.count.timer_print_progress = 0;
        if (pop.world_ID == 0) {
            IOFormat fmt(19, 0, "  ", "  ", "", "", "", "");
            cout << fixed << setprecision(19);
            if (pop.obj_fun.id >= 0) { cout << "ID: " << pop.obj_fun.id << " " << pop.obj_fun.name << " "; }
            cout << "Generation... " << setw(7) << pop.count.generation
                 << " | Current Best: " << setw(22) << champion_fitness
                 << " | diff: "         << setw(22) << latest_history_diff
//                 << endl << champion_parameters.transpose()
                 << endl << flush;
        }
    }
}

void hall_of_fame::print_final_result(population &pop){

    if (pop.world_ID == 0){
        if (pop.obj_fun.id >= 0){cout << "ID: " << pop.obj_fun.id << " "<< pop.obj_fun.name << " ";}
        cout << "Best Parameters give H = : "
             << champion_fitness
             << endl;
    }

}

void hall_of_fame::store_champion(population &pop) {
    if (pop.count.timer_store_best++ > rate_store_best) {
        pop.count.timer_store_best = 0;
        update_global_champion(pop);
        previous_champions_fitness.push_back(champion_fitness);
    }
}

bool hall_of_fame::converged(population &pop) {
    if (pop.count.timer_check_conv++ > rate_check_conv) {
        pop.count.timer_check_conv = 0;
        //How many generations do we have stored?
        int gen = rate_store_best * (int)previous_champions_fitness.size();
        //How many generations back do we want to look? num_check_history.
        //From which generation should we check?  gen - num_check_history
        //Which index does this correspond to?
        int from =  (gen - num_check_history) / rate_store_best;
        int n = (int)previous_champions_fitness.size() - from - 1 ;
        if (from > 0 && n > 0){
            ArrayXd fitness_history = Map<ArrayXd>(previous_champions_fitness.data(), previous_champions_fitness.size());
            latest_history_diff = (fitness_history.segment(from-1, n) - fitness_history.segment(from,n)).mean();
//            cout << fitness_history.segment(from, n).transpose() << endl;
            pop.current_diff = latest_history_diff;
        }
        return fabs(latest_history_diff) < pop.obj_fun.tolerance;
    }else{
        return false;
    }
}

//        if(fitness_history.size() <= count.store_counter) {
//            fitness_history.conservativeResize(2 * max(1,count.store_counter));
//        }
//        double champ = champion_fitness();
//        if (count.store_counter > 0){
//            if(champ == fitness_history(count.store_counter-1)){
//                count.store_last_since++;
//            }else{
//                fitness_history(count.store_counter) = champ;
//                count.store_counter++;
//                count.store_last_since = 0;
//            }
//            if(count.store_last_since == 5){ //Override
//                fitness_history(count.store_counter) = champ;
//                count.store_counter++;
//                count.store_last_since = 0;
//            }
//        }else{
//            fitness_history(count.store_counter) = champ;
//            count.store_counter++;
//            count.store_last_since = 0;
//        }
//    }
//
//}


