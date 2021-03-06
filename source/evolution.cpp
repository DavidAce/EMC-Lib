
#include "evolution.hpp"

void evolve(population &pop) {
	//This function selects either mutation type operators or crossover type.
	//Then does an exchange operator, and finally finds the best guys in the population
	pop.generation++;
	int dice;
	//Select mutation or crossover
	if (uniform_double_1() < qm) {
        if (uniform_double_1() < qma) {
			mutation_elite(pop);
		}
		else {
			mutation(pop);
		}
	}
	else {
		dice = uniform_integer(0, 2);
		if 		(dice == 0)	 {if(uniform_double_1() < qe){crossover_elite(pop); }else{crossover(pop);}}
		else if (dice == 1)	 { crossover_smartCopy(pop); }
		else				 {crossover_snooker  (pop); }


	}
    exchange(pop);
    find_elite(pop);

}


void find_lowest_guy(const vector<personality> &guys, long double &lowest_H, unsigned int &lowest_i ){
    lowest_H = guys[0].H;
    lowest_i = 0;
    for (unsigned int i = 0; i < guys.size(); i++){
        if (guys[i].H < lowest_H){
            lowest_H = guys[i].H;
            lowest_i = i;
        }
    }
}

void find_lowest_guy_excluding(const vector<personality> &guys, long double &lowest_H, unsigned int &lowest_i, std::set<unsigned int> &excluded ){
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

void exchange(population &pop) {
    int i, j, ex; //ex is the exchange iteration
    long double re, dH, dt_inv;
    for (ex = 0; ex < N; ex++) {
        i = uniform_integer( 0, N - 1);										//Chose a guy...
        j = i == 0 ? 1 : (i== N-1? N-2: i + 2 * uniform_integer( 0, 1) - 1);	//And a neighbor
        dH = pop.guys[i].H - pop.guys[j].H;
        dt_inv = 1 / pop.guys[i].t - 1 / pop.guys[j].t;
        re = exp(dH*dt_inv);

        if (uniform_double_1() < fmin(1, re)) {
            //Save the date of birth
            pop.newguys[j].born = pop.guys[i].born;
            pop.newguys[i].born = pop.guys[j].born;
            //Swap the guys on the ladder...
            pop.copy(pop.guys[i], pop.newguys[j]);
            pop.copy(pop.guys[j], pop.newguys[i]);
            //But keep the temperature of position i as t_i,and j as t_j
            pop.guys[i].t = pop.newguys[i].t;
            pop.guys[j].t = pop.newguys[j].t;
            //Now the list of newguys should be updated to be identical to guys
            pop.copy(pop.newguys[i], pop.guys[i]);
            pop.copy(pop.newguys[j], pop.guys[j]);
        }
    }
}
void migration(species &sp) {
    //First select randomly m out of M pobulations, subject to migration
    //Then select randomly  n (n=1?) out of N guys to migrate
    //Accept the new population based on the new population fitness score
    int n;										//Guy chosen to migrate
    long double dH;									//Probability of migration
    int m = uniform_integer(2, M);		//How many populations to migrate
    ArrayXi s_pop(m);							//Array of sender populations
    ArrayXi r_pop(m);							//Array of receiving populations
    rndChoice(s_pop.data(), m, M);				//Fill the array of migrators with m random populations

    for (int i = 0; i < m; i++) {
        r_pop(i) = i < m - 1 ? s_pop(i + 1) : s_pop(0) ;						//Define receiving populations cyclically
        n = uniform_integer( 0, N-1);									//Who will be the lucky one?
        sp.copy(sp.pop[r_pop(i)].newguys[n], sp.pop[s_pop(i)].newguys[n]);	//Inject a guy into another population
        dH = sp.pop[r_pop(i)].newguys[n].H - sp.pop[r_pop(i)].guys[n].H;
        if (dH < 0 || exp(-dH / sp.pop[r_pop(i)].newguys[n].t) > uniform_double_1()) {
            sp.pop[r_pop(i)].newguys[n].born = sp.count.generation;
            sp.copy(sp.pop[r_pop(i)].guys[n], sp.pop[r_pop(i)].newguys[n]);
        }
        else {
            //Revert changes in newguys,  i.e sync them for the next round
            sp.copy(sp.pop[r_pop(i)].newguys[n], sp.pop[r_pop(i)].guys[n]);
        }
    }

}
void insertguy(population &pop, int from, int to) {
    //Push out the worst bestguy and move everybody up until "to"
    for (int i = 0; i < to; i++) {
        pop.copy(pop.bestguys[i], pop.bestguys[i + 1]);
    }
    //Insert the exceptional guy into the illustrious group of best guys
    pop.copy(pop.bestguys[to], pop.guys[from]);
}

void find_elite(population &pop) {
    //Here we attempt to find someone better in guys that is better than any guy in best_guys
    std::set<unsigned int> tried_guys;
    long double lowest_H;
    unsigned int lowest_i;
    int j = 0; //Counts how many have been tried
    while (j < N_best) {
        lowest_H = pop.guys[0].H;
		lowest_i = 0;
        //Find the best guy yet among guys
        find_lowest_guy_excluding(pop.guys, lowest_H,lowest_i, tried_guys);
        tried_guys.insert(lowest_i);
        j++;
        //By now we should have a winner, check if he's better than any elite
        //We have already tried to match j guys, so we can start from j?
        for (unsigned int i = N_best ;  i-- > 0;) {
            if (lowest_H < pop.bestguys[i].H) {
                insertguy(pop, lowest_i, i);
                break;
            }
        }

    }
}


//void roulette_select(vector<personality> &guys, ArrayXi &selected, long double &Z, double s) {
//    //Selected 0 is the guy with good fitness (lowest), selected 1 is random
//
//    //double total_H = 0;
//    ArrayXd roulette(N);
//    double lucky_number;
//    //Make a roulette wheel
//    Z = 0;
//    for (int i = 0; i < N; i++) {
//        Z += static_cast<long double>(exp(-guys[i].H / s));
//        roulette(i) = static_cast <double> (Z); //Cumulative sum - Low fitness gives large area on roulette
//    }
//    lucky_number = uniform_double( 0.0 , static_cast<double>(Z));
//    for (int i = 0; i < N; i++) {
//        if (lucky_number <= roulette(i)) {
//            selected(0) = i;
//            break;
//        }
//    }
//    selected(1) = selected(0);
//    while (selected(1) == selected(0)) {
//        selected(1) = uniform_integer( 0, N - 1);
//    }
//}

void roulette_select(const vector<personality> &guys, Array2i &selected, long double &Z, long double lowest_H) {
	//Selected 0 is the guy with good fitness (lowest), selected 1 is random

	Array<long double, N, 1> roulette;

	long double lucky_number;
	//Make a roulette wheel
    //Start by finding the best guy in the list
	Z = 0;
	for (int i = 0; i < N; i++) {
		Z += exp(-(guys[i].H - lowest_H));
		roulette(i) = Z; //Cumulative sum - Low fitness gives large area on roulette
//        allZ(i) = exp(-(guys[i].H - lowest_H));
//        allH(i) = (guys[i].H - lowest_H);
//        allt(i) = (guys[i].t);
	}
//    #pragma omp critical
//    {
//        cout << lowest_H << endl;
//        cout << allH.transpose() << endl;
//        cout << allt.transpose() << endl;
//        cout << allZ.transpose() << endl << endl;
//    }
	lucky_number = uniform_double( (long double) 0.0 , Z);
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


inline void bitselector_smartCopy(population &pop, Array2i &selected, Array4i &n_ab_XY, Array4i &n_ab_YX) {
	//n_ab_XY(0): n_11 = guys: SAME | newguys: SAME
	//n_ab_XY(1): n_12 = guys: SAME | newguys: DIFF
	//n_ab_XY(2): n_21 = guys: DIFF | newguys: SAME
	//n_ab_XY(3): n_22 = guys: DIFF | newguys: DIFF

	//n_ab_YX(0): n_11 = newguys: SAME | guys: SAME
	//n_ab_YX(1): n_12 = newguys: SAME | guys: DIFF
	//n_ab_YX(2): n_21 = newguys: DIFF | guys: SAME
	//n_ab_YX(3): n_22 = newguys: DIFF | guys: DIFF
	
	//#pragma omp parallel for 
	for (int i = 0; i < genomeLength; i++) {
		if (pop.guys[selected(0)].genome(i) == pop.guys[selected(1)].genome(i)) {
			//If parents loci are the same: Copy the values
			pop.newguys[selected(0)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
			pop.newguys[selected(1)].genome.copy_loci(i, pop.guys[selected(1)].genome(i));

			//Independently reverse with probability p0
			for (int j = 0; j < 2; j++) {
				if (uniform_double_1() < P0) {
					pop.newguys[selected[j]].genome.flip_loci(i);
				}
			}

			if (pop.newguys[selected(0)].genome(i) == pop.newguys[selected(1)].genome(i)) {
				//#pragma omp atomic
				n_ab_XY(0)++; //Offsprings have identical bit
				//#pragma omp atomic
				n_ab_YX(0)++; //Parents have identical bit
			}
			else {
				//#pragma omp atomic
				n_ab_XY(1)++; //Offsprings have different bit
				//#pragma omp atomic
				n_ab_YX(2)++; //Parents have different bit
			}

		}
		else {
			//If parents loci are different: Copy the values
			pop.newguys[selected(0)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
			pop.newguys[selected(1)].genome.copy_loci(i, pop.guys[selected(1)].genome(i));
			
			//Reverse with probabilities p1 and p2 respectively
			if (uniform_double_1() < P1) {
				pop.newguys[selected(0)].genome.flip_loci(i);
			}
			if (uniform_double_1() < P2) {
				pop.newguys[selected(1)].genome.flip_loci(i);
			}

			if (pop.newguys[selected(0)].genome(i) == pop.newguys[selected(1)].genome(i)) {
				//#pragma omp atomic
				n_ab_XY(2)++; //Offsprings have identical bit
				//#pragma omp atomic
				n_ab_YX(1)++; //Offsprings have identical bit
			}
			else {
				//#pragma omp atomic
				n_ab_XY(3)++; //Offsprings have different bit
				//#pragma omp atomic
				n_ab_YX(3)++; //Parents have different bit

			}

		}
	}
}

void mutation(population &pop) {
	int mutantGenes;	//Number of points to mutate
	int mutant;			//Which guy to mutate
	long double dH;
	//#pragma omp parallel for private(mutantGenes, mutant,dH)
	for (int i = 0; i < N; i++) {
		mutantGenes =  uniform_integer( 1, genomeLength - 1);
		ArrayXi loci(mutantGenes);
		mutant = i;// uniform_integer( 0, N - 1);
		rndChoice(loci.data(), mutantGenes, genomeLength);				//Choose locus to mutate
		pop.newguys[mutant].genome.flip_loci(loci);	                    //Flip bits
		pop.getFitness(pop.newguys[mutant]);							//Get Fitness score
		
		//Perform Metropolis
		dH = pop.newguys[mutant].H - pop.guys[mutant].H;		//used to decide if we accept the new guy or not.
		if (dH < 0 || exp(-dH / pop.newguys[mutant].t) > uniform_double_1()) {
			pop.newguys[mutant].born = pop.generation;
			pop.copy(pop.guys[mutant], pop.newguys[mutant]);

		}
		else { 	
			//Revert changes in newguys,  i.e sync them for the next round
			pop.copy(pop.newguys[mutant], pop.guys[mutant]);
		}
	}


}

void mutation_elite(population &pop) {
	//int mutantGenes;	//Number of points to mutate
	int mutant;			//which guy to mutate
	int elite_mutant;	//Which elite guy to receive wisdom from
	long double dH;
	//#pragma omp parallel for private( mutant, elite_mutant,dH)
	for (int i = 0; i < N; i++) {
		//mutantGenes = 1;// uniform_integer( 1, genomeLength - 1);
		int loci = uniform_integer( 1, genomeLength - 1);
		mutant = i; //uniform_integer( 0, N - 1);
		//Fill loci with mutantGenes genome points to be mutated
		//Copy an elite guy to a new guy to be a guinnea pig
		elite_mutant = uniform_integer( 0, N_best - 1);
		pop.copy(pop.newguys[mutant].genome, pop.bestguys[elite_mutant].genome);	//Copy DNA only!
		pop.newguys[mutant].genome.flip_loci(loci);									//Flip bits	
		pop.getFitness(pop.newguys[mutant]);    									//Get Fitness score
		
		//Perform metropolis
		dH = pop.newguys[mutant].H - pop.guys[mutant].H;		//used to decide if we accept the new guy or not.
		if (dH < 0 || exp(-dH / pop.newguys[mutant].t) > uniform_double_1()) {
			pop.newguys[mutant].born = pop.generation;
			pop.copy(pop.guys[mutant], pop.newguys[mutant]);
		}
		else {
			//Revert changes in newguys, i.e sync them for the next round
			pop.copy(pop.newguys[mutant], pop.guys[mutant]);
		}

	}
}

void crossover(population &pop) {
	int nMatings = (int)(0.2*N);
	Array2i selected(2);
    Array<long double, 2,1> expX;
    Array<long double, 2,1> expY;
//	ArrayXd expX(2);
//	ArrayXd expY(2);
	int crossoverPoint;
	long double rc;
    long double dHt0, dHt1;
	long double PXX, PYY;			//Selection probabilities
	long double PXY = 1, PYX = 1;	//Generating probabilities (PXY = PYX for this operator)
	long double TXY, TYX;			//Transition probability
    long double lowest_H;
    unsigned int lowest_i;
	long double ZX, ZY;		//Sum of Boltzmann-weights for current (X) and offspring(Y) populations

	for (int matings = 0; matings < nMatings; matings++) {
        ZX=0;
        ZY=0;
        find_lowest_guy(pop.guys,lowest_H,lowest_i);
		roulette_select(pop.guys, selected, ZX, lowest_H); //Selected 0 will be good, selected 1 random
        expX(0) = exp(-(pop.guys[selected(0)].H-lowest_H)); //good guy
		expX(1) = exp(-(pop.guys[selected(1)].H-lowest_H)); //bad guy
		PXX = ( 1 / ((N - 1)*ZX)*(expX(0) + expX(1))); //P((xi,xj) | x)
		
		//Now mate the newguys to create offsprings
		crossoverPoint = uniform_integer( 1, genomeLength-1);
		//cout << selected.transpose() << " " << pop.newguys[selected(0)].H << endl;
		for (int i = crossoverPoint; i < genomeLength; i++) {
			pop.newguys[selected(0)].genome.copy_loci(i, pop.guys[selected(1)].genome(i));
			pop.newguys[selected(1)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
		}
		for (int i = 0; i < 2; i++) {
			pop.getFitness(pop.newguys[selected(i)]);
		}
        find_lowest_guy(pop.newguys,lowest_H,lowest_i);
        expY(0) = exp(-(pop.newguys[selected(0)].H-lowest_H));
		expY(1) = exp(-(pop.newguys[selected(1)].H-lowest_H));
		for (int i = 0; i < N; i++) {
			ZY += exp(-(pop.newguys[i].H - lowest_H));
		}
		PYY = (1 / ((N - 1)*ZY)*(expY(0) + expY(1))); //P((yi,yj) | y) selection probability
		dHt0 = (pop.newguys[selected(0)].H - pop.guys[selected(0)].H) / pop.guys[selected(0)].t; //good 
		dHt1 = (pop.newguys[selected(1)].H - pop.guys[selected(1)].H) / pop.guys[selected(1)].t; //bad 

		TXY = PXX*PXY;
		TYX = PYY*PYX;

		rc = exp(-dHt0 - dHt1)*TXY / TYX;
//        if (dHt0 < 10 || dHt1 < 10){
//            cout << "rc = " <<rc << " dHt0 = " << dHt0 << " dHt1 = " << dHt1 << " "  << TXY << " " << TYX <<  endl;
//            cout << selected.transpose() << endl;
//            cout << "  0 New " <<pop.newguys[0].genome.parameters.transpose();
//            cout << "    Old " <<pop.guys[0].genome.parameters.transpose() << " t: "<<pop.guys[selected(0)].t<<  endl;
//            cout << "  1 New " <<pop.newguys[1].genome.parameters.transpose();
//            cout << "    Old " <<pop.guys[1].genome.parameters.transpose() << " t: "<<pop.guys[selected(1)].t <<  endl;
//        }
		//cout << "Selected:	" << selected.transpose()<< "	P: " << rc << " dHt0: " << dHt0 << " dHt1: " << dHt1 << endl;
		//Accept or reject
		if (uniform_double_1() < fmin(1, rc)) {
			pop.newguys[selected(0)].born = pop.generation;
			pop.newguys[selected(1)].born = pop.generation;
			pop.copy(pop.guys[selected(0)], pop.newguys[selected(0)]);
			pop.copy(pop.guys[selected(1)], pop.newguys[selected(1)]);
		}
		else {
			//Revert changes on newguys
			pop.copy(pop.newguys[selected(0)], pop.guys[selected(0)]);
			pop.copy(pop.newguys[selected(1)], pop.guys[selected(1)]);
		}
		//Make sure newguys are up to speed on newest events

		

		
	/*	cout << "Parameters 0:	" << pop.guys[selected(0)].genome.parameters.transpose() << endl;
		cout << "Parameters 1:	" << pop.guys[selected(1)].genome.parameters.transpose() << endl;
		getchar();*/
	}
}



void crossover_elite(population &pop) {
	//Start roulette selection
	int matings;
	int nMatings = (int)(0.2*N);
	Array2i selected(2);
	int crossoverPoint;
    long double ZX,ZY;
	long double rc, dHt0, dHt1;
	long double PXX, PYY;			//Selection probabilities
	long double PXY, PYX;	        //Generating probabilities (PXY = PYX for this operator)
	long double TXY, TYX;			//Transition probability
    Array<long double, 2,1> expX;
    Array<long double, 2,1> expY;
    long double lowest_H;
    unsigned int lowest_i;
	int random_bestguy; //Guy to inject into selected[1]
	for (matings = 0; matings < nMatings; matings++) {
        ZX = 0;
        ZY = 0;
        //Selected 0 is a boltzmann "just good" guy. Selected 1 is a random any guy.
        find_lowest_guy(pop.guys,lowest_H,lowest_i);
        roulette_select(pop.guys, selected, ZX, lowest_H); //Selected 0 will be good, selected 1 random

        //Let the newguy selected[1] impersonate a random bestguy but keep the temperature. Then newguy gets superpowers
		random_bestguy = uniform_integer( 0, N_best-1);
		pop.copy(pop.newguys[selected(1)], pop.bestguys[random_bestguy]);
		pop.newguys[selected(1)].born = pop.guys[selected(1)].born;		//Keep date of birth count
		pop.newguys[selected(1)].t = pop.guys[selected(1)].t;			//Keep temperature

		//Now selected(1) is some guy high on the temperature ladder with amazing bestguy-genes and fitness
		expX(0) = exp(-(pop.newguys[selected(0)].H - lowest_H)); //good guy
		expX(1) = exp(-(pop.newguys[selected(1)].H - lowest_H)); //gooder guy
        PXX = ( 1 / ((N - 1)*ZX)*(expX(0) + expX(1))); //P((xi,xj) | x)

        //Make sure the guy selected(0) and random_bestguys aren't already the same guy, or else they will have identical
        //Offspring. Basically we would be copying back a bestguy into our list.
        if ( pop.guys[selected(0)]. H == pop.bestguys[random_bestguy].H){
            //Revert changes
            pop.copy(pop.newguys[selected[0]], pop.guys[selected[0]]);
            pop.copy(pop.newguys[selected[1]], pop.guys[selected[1]]);
            continue;
        }
		//Now mate to create offspring. Let a bestGuy inject DNA in this process!
		crossoverPoint = uniform_integer( 1, genomeLength - 1);
		for (int i = 0; i < genomeLength; i++) {
			if (i < crossoverPoint) {
				pop.newguys[selected(0)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
				pop.newguys[selected(1)].genome.copy_loci(i, pop.bestguys[random_bestguy].genome(i));
			}
			else {
				pop.newguys[selected(0)].genome.copy_loci(i, pop.bestguys[random_bestguy].genome(i));
				pop.newguys[selected(1)].genome.copy_loci(i, pop.guys[selected(0)].genome(i));
			}
		}

		//From now on the newguys are offsprings. The parents are a good guy selected 1 and a bestguy random_bestguy. Adoptive father is guy selected 1
		
		for (int i = 0; i < 2; i++) {
			pop.getFitness(pop.newguys[selected(i)]);
		}
        find_lowest_guy(pop.newguys,lowest_H,lowest_i);
		expY(0) = exp(-(pop.newguys[selected(0)].H - lowest_H));
		expY(1) = exp(-(pop.newguys[selected(1)].H - lowest_H));
		for (int i = 0; i < N; i++) {
			ZY += exp(-(pop.newguys[i].H - lowest_H));
		}
		PYY = (1 / ((N - 1)*ZY)*(expY(0) + expY(1))); //P((xi,xj) | x)
        PXY = 1;
        PYX = PXY;
		TXY = PXX*PXY;
		TYX = PYY*PYX;

		dHt0 = (pop.newguys[selected(0)].H - pop.guys[selected(0)].H) / pop.guys[selected(0)].t; //good 
		dHt1 = (pop.newguys[selected(1)].H - pop.guys[selected(1)].H) / pop.guys[selected(1)].t; //bad 


		rc = exp(-dHt0 - dHt1)*TXY / TYX;
//        if (true){
//            cout << "rc = " <<rc << " dHt0 = " << dHt0 << " dHt1 = " << dHt1 << " "  << TXY << " " << TYX <<  endl;
//            cout << selected.transpose() << endl;
//            cout << "  0 New " <<pop.newguys[0].genome.parameters.transpose();
//            cout << "    Old " <<pop.guys[0].genome.parameters.transpose() << " t: "<<pop.guys[selected(0)].t<<  endl;
//            cout << "  1 New " <<pop.newguys[1].genome.parameters.transpose();
//            cout << "    Old " <<pop.guys[1].genome.parameters.transpose() << " t: "<<pop.guys[selected(1)].t <<  endl;
//        }

		//Accept or reject
		if (uniform_double_1() < fmin(1, rc)) {
			//Keep children.
			pop.newguys[selected(0)].born = pop.generation;
			pop.newguys[selected(1)].born = pop.generation;
			pop.copy(pop.guys[selected[0]], pop.newguys[selected[0]]);
			pop.copy(pop.guys[selected[1]], pop.newguys[selected[1]]);
		}
		else {
			//Keep parents.
			pop.copy(pop.newguys[selected[0]], pop.guys[selected[0]]);
			pop.copy(pop.newguys[selected[1]], pop.guys[selected[1]]);
		}

	}

}

void crossover_smartCopy(population &pop) {
	//Start roulette selection																							
	int matings;
	int nMatings = (int)(0.2*N);
	Array2i selected;
    long double ZX, ZY;
    long double rc, dHt0, dHt1;
	long double PXX, PYY; //Selection probabilities
	long double PXY, PYX; //Generating probabilities (PXY != PYX for this operator)
	long double TXY, TYX; //Transition probability
    Array<long double, 2,1> expX;
    Array<long double, 2,1> expY;
    long double lowest_H;
    unsigned int lowest_i;
	Array4i n_ab_XY; //Exponents for smartCopy probabilities. Note that n_ab.sum() = genomeLength
	Array4i n_ab_YX; //Exponents for smartCopy probabilities. Note that n_ab.sum() = genomeLength
	for (matings = 0; matings < nMatings; matings++) {
        ZX=0;
        ZY=0;
        n_ab_XY.fill(0);
		n_ab_YX.fill(0);
        find_lowest_guy(pop.guys,lowest_H, lowest_i);
        roulette_select(pop.guys, selected, ZX, lowest_H); //Selected 0 will be good, selected 1 random
        expX(0) = exp(-(pop.guys[selected(0)].H-lowest_H)); //good guy
        expX(1) = exp(-(pop.guys[selected(1)].H-lowest_H)); //bad guy
        PXX = ( 1 / ((N - 1)*ZX)*(expX(0) + expX(1))); //P((xi,xj) | x)

        //Now mate the newguys to create offsprings
		bitselector_smartCopy(pop, selected, n_ab_XY,n_ab_YX);
		
		//Now we need generating probabilities PXY and PYX
		PXY = 1;
		PYX = 1;
		for (int i = 0; i < 4; i++) {
			PXY *= pow(p_matrix[i], n_ab_XY(i)/genomeLength); //Convert exponent to frequency
			PYX *= pow(p_matrix[i], n_ab_YX(i)/genomeLength); //Convert exponent to frequency
		}
		for (int i = 0; i < 2; i++) {
			pop.getFitness(pop.newguys[selected(i)]);
		}
        find_lowest_guy(pop.newguys,lowest_H,lowest_i);
        expY(0) = exp(-(pop.newguys[selected(0)].H-lowest_H));
        expY(1) = exp(-(pop.newguys[selected(1)].H-lowest_H));
        for (int i = 0; i < N; i++) {
            ZY += exp(-(pop.newguys[i].H - lowest_H));
        }
		PYY = (1 / ((N - 1)*ZY)*(expY(0) + expY(1))); //P((yi,yj) | y) selection probability
	    dHt0 = (pop.newguys[selected(0)].H - pop.guys[selected(0)].H) / pop.guys[selected(0)].t; //good
		dHt1 = (pop.newguys[selected(1)].H - pop.guys[selected(1)].H) / pop.guys[selected(1)].t; //bad 
		
		TXY = PXX*PXY;
		TYX = PYY*PYX;

		rc = exp(-dHt0 - dHt1)*TXY / TYX;
//        cout << rc << " dHt0 = " << dHt0 << " dHt1 = " << dHt1 << " " << TXY << " " << TYX <<  endl;
		//Accept or reject
		if (uniform_double_1() < fmin(1, rc)) {
			//Keep Children
			pop.newguys[selected(0)].born = pop.generation;
			pop.newguys[selected(1)].born = pop.generation;
			pop.copy(pop.guys[selected[0]], pop.newguys[selected[0]]);
			pop.copy(pop.guys[selected[1]], pop.newguys[selected[1]]);

		}
		else {
			//Keep parents
			pop.copy(pop.newguys[selected[0]], pop.guys[selected[0]]);
			pop.copy(pop.newguys[selected[1]], pop.guys[selected[1]]);
		}

	}
}

void crossover_snooker(population &pop) {
	//Perhaps the distribution f(r) can be done with a roulette?
	//f(r) = exp(H(x+re))/integral_-r_min^r_max exp(H(x+re)) typ?
	//Notation from snooker crossover
	//Choose 3 guys from pop
	
	Array2i selected; //Guy 0 is "x_i", 1 is x_j, "anchor";
	rndChoice(selected.data(), 2, N - 1);
	// ArrayXd &p0 = pop.guys[selected(0)].genome.parameters; //Reference to guy 0, 
	// ArrayXd &p1 = pop.guys[selected(1)].genome.parameters; //Reference to guy 1, who is lower on the ladder (better)
	pop.line.Through(pop.guys[selected(0)].genome.parameters, pop.guys[selected(1)].genome.parameters);
	//Distance in terms of r between guys
	//double distance = pop.line.distance(pop.guys[selected(1)].genome.parameters,pop.guys[selected(0)].genome.parameters);
	//Vary r to find the walls of the parameter domain
    long double r_max = pop.line.line_max();
    long double r_min = pop.line.line_min();
	Array<long double, Dynamic, 1> r_point(r_num);
    if( isinf(r_min) || isinf(r_max)){
        return;
    }
    if( r_min * r_max > 0){
        return;
    }
	for (int i = 0; i < r_num; i++) {
		//r_point(i) = uniform_double( fmin(r_min,r_max), fmax(r_min,r_max));
		r_point(i) = gaussian_truncated(	//Create gaussians points centered around the good performer
										fminl(r_min, r_max),
										fmaxl(r_min, r_max),
										1.0,
										0.3);
		pop.snookerGuys[i].genome.set_parameters(pop.line.pointAt(r_point(i)));
		pop.getFitness(pop.line.pointAt(r_point(i)), pop.snookerGuys[i]);
		pop.snookerGuys[i].t = pop.guys[selected(0)].t;
	}
	//Time to make a roulette to see which "r" is chosenpop.

	Array<long double, r_num,1> roulette;
	double lucky_number = uniform_double_1();
    long double lowest_H;
    unsigned int lowest_i;

    find_lowest_guy(pop.snookerGuys, lowest_H,lowest_i);
	//Make a roulette wheel
	long double Z = 0;

    for (int i = 0; i < r_num; i++) {
		Z += exp(-(pop.snookerGuys[i].H - lowest_H)); //Area on roulette wheel proportional to Boltzmann weights
		roulette(i) = Z; //Cumulative sum - Low fitness gives large area on roulette
	}
	roulette /= Z; //Normalize to a proper distribution
	for (int i = 0; i < r_num; i++) {
		if (lucky_number <= roulette(i)) {
			pop.guys[selected(0)].genome.set_parameters(pop.snookerGuys[i].genome.parameters); //Copy his shit
			pop.guys[selected(0)].H 	= pop.snookerGuys[i].H;
			pop.guys[selected(0)].born 	= pop.generation;
			break;
		}
	}
}


