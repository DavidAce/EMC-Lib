
#include "evolution.hpp"

void evolve(population &pop) {
	//This function selects either mutation type operators or crossover type.
	//Then does an exchange operator, and finally finds the best guys in the population
	int dice;
	//Select mutation or crossover
	if (uniform_double_1() < qm) {
		dice = uniform_integer(0,1);
		if 		(dice == 0)	 {if(uniform_double_1() < qma){mutation_elite(pop); }else{mutation(pop);}}
		else   				 { mutation_realspace(pop); }

	}
	else {
		dice = uniform_integer(1, 2);
		if 		(dice == 1)	 {if(uniform_double_1() < qe){crossover_elite(pop); }else{crossover(pop);}}
		else if (dice == 0)	 { crossover_smartCopy(pop); }
		else				 { crossover_snooker  (pop); }


	}
    exchange(pop);
    pop.find_elite();
    pop.count.generation++;
}



void exchange(population &pop) {
    int i, j, ex; //ex is the exchange iteration
    double re, dH, dt_inv;
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






inline void bitselector_smartCopy(population &pop, Array4i &n_ab_XY, Array4i &n_ab_YX) {
	//n_ab_XY(0): n_11 = guys: SAME | newguys: SAME
	//n_ab_XY(1): n_12 = guys: SAME | newguys: DIFF
	//n_ab_XY(2): n_21 = guys: DIFF | newguys: SAME
	//n_ab_XY(3): n_22 = guys: DIFF | newguys: DIFF

	//n_ab_YX(0): n_11 = newguys: SAME | guys: SAME
	//n_ab_YX(1): n_12 = newguys: SAME | guys: DIFF
	//n_ab_YX(2): n_21 = newguys: DIFF | guys: SAME
	//n_ab_YX(3): n_22 = newguys: DIFF | guys: DIFF
	
	for (int i = 0; i < genomeLength; i++) {
		if (pop.guys[pop.selected(0)].genome(i) == pop.guys[pop.selected(1)].genome(i)) {
			//If parents loci are the same: Copy the values
			pop.newguys[pop.selected(0)].genome.copy_loci(i, pop.guys[pop.selected(0)].genome(i));
			pop.newguys[pop.selected(1)].genome.copy_loci(i, pop.guys[pop.selected(1)].genome(i));

			//Independently reverse with probability p0
			for (int j = 0; j < 2; j++) {
				if (uniform_double_1() < P0) {
					pop.newguys[pop.selected[j]].genome.flip_loci(i);
				}
			}

			if (pop.newguys[pop.selected(0)].genome(i) == pop.newguys[pop.selected(1)].genome(i)) {
				n_ab_XY(0)++; //Offsprings have identical bit
				n_ab_YX(0)++; //Parents have identical bit
			}
			else {
				n_ab_XY(1)++; //Offsprings have different bit
				n_ab_YX(2)++; //Parents have different bit
			}

		}
		else {
			//If parents loci are different: Copy the values
			pop.newguys[pop.selected(0)].genome.copy_loci(i, pop.guys[pop.selected(0)].genome(i));
			pop.newguys[pop.selected(1)].genome.copy_loci(i, pop.guys[pop.selected(1)].genome(i));
			
			//Reverse with probabilities p1 and p2 respectively
			if (uniform_double_1() < P1) {
				pop.newguys[pop.selected(0)].genome.flip_loci(i);
			}
			if (uniform_double_1() < P2) {
				pop.newguys[pop.selected(1)].genome.flip_loci(i);
			}

			if (pop.newguys[pop.selected(0)].genome(i) == pop.newguys[pop.selected(1)].genome(i)) {
				n_ab_XY(2)++; //Offsprings have identical bit
				n_ab_YX(1)++; //Offsprings have identical bit
			}
			else {
				n_ab_XY(3)++; //Offsprings have different bit
				n_ab_YX(3)++; //Parents have different bit

			}

		}
	}
}

//void mutation(population &pop) {
//	int mutant;			//Which guy to mutate
//	int mutantGenes =  uniform_integer( 2, genomeLength - 1);//Number of points to mutate
//	ArrayXi loci(mutantGenes);
//	double dH;
//	for (int i = 0; i < N; i++) {
//
//		mutant = i;// uniform_integer( 0, N - 1);
//		rndChoice(loci.data(), mutantGenes, genomeLength);				//Choose locus to mutate
//		pop.newguys[mutant].genome.flip_loci(loci);	                    //Flip bits
//		pop.getFitness(pop.newguys[mutant]);							//Get Fitness score
//
//		//Perform Metropolis
//		dH = pop.newguys[mutant].H - pop.guys[mutant].H;		//used to decide if we accept the new guy or not.
//		if (dH < 0 || exp(-dH / pop.newguys[mutant].t) > uniform_double_1()) {
//			pop.newguys[mutant].born = pop.count.generation;
//			pop.copy(pop.guys[mutant], pop.newguys[mutant]);
//
//		}
//		else {
//			//Revert changes in newguys,  i.e sync them for the next round
//			pop.copy(pop.newguys[mutant], pop.guys[mutant]);
//		}
//	}
//
//
//}
//
//void mutation_elite(population &pop) {
//	int mutantGenes;	//Number of points to mutate
//	int mutant;			//which guy to mutate
//	int elite_mutant;	//Which elite guy to receive wisdom from
//
//	double dH; 			//Fitness difference
//	ArrayXi loci(genomeLength);
////	ArrayXi sub_loci(geneLength);
//	for (int i = 0; i < N; i++) {
//		mutant = i;
//		elite_mutant = uniform_integer( 0, N_best - 1);
//		pop.copy(pop.newguys[mutant].genome, pop.bestguys[elite_mutant].genome);	//Copy DNA only!
////		cout <<endl<< pop.newguys[mutant].genome << endl << endl;
//
//		mutantGenes =  uniform_integer( 0, genomeLength - 1);
//		rndChoice(loci.data(), mutantGenes, genomeLength);
//		pop.newguys[mutant].genome.flip_loci(loci.head(mutantGenes));
//		pop.getFitness(pop.newguys[mutant]);    									//Get Fitness score
////		cout << pop.newguys[mutant].genome << endl << endl <<  endl;
//
//		//Perform metropolis
//		dH = pop.newguys[mutant].H - pop.guys[mutant].H;		//used to decide if we accept the new guy or not.
//		if (dH < 0 || exp(-dH / pop.newguys[mutant].t) > uniform_double_1()) {
//			pop.newguys[mutant].born = pop.count.generation;
//			pop.copy(pop.guys[mutant], pop.newguys[mutant]);
//		}
//		else {
//			//Revert changes in newguys, i.e sync them for the next round
//			pop.copy(pop.newguys[mutant], pop.guys[mutant]);
//		}
//	}
//}
void mutation(population &pop) {
	int mutantGenes;	//Number of points to mutate
	int mutant;			//Which guy to mutate
	double dH;
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
			pop.newguys[mutant].born = pop.count.generation;
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
	double dH;
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
			pop.newguys[mutant].born = pop.count.generation;
			pop.copy(pop.guys[mutant], pop.newguys[mutant]);
		}
		else {
			//Revert changes in newguys, i.e sync them for the next round
			pop.copy(pop.newguys[mutant], pop.guys[mutant]);
		}
	}
}


void mutation_realspace(population &pop) {
	int mutantGenes;	//Number of genes to mutate
	int mutant;			//which guy to mutate
	double dH; 			//Fitness difference
	for (int i = 0; i < N; i++) {
		mutant = i;
		mutantGenes = uniform_integer(1,nGenes);
		for (int j = 0; j < mutantGenes ; j++){
			int k = uniform_integer( 0, nGenes - 1);
            double s = pop.current_diff*fabs(fabs(pop.obj_fun.upper_bound(k) - pop.obj_fun.lower_bound(k)));
            pop.newguys[mutant].genome.parameters(k) = EMC_rnd::gaussian_truncated(pop.obj_fun.lower_bound(k),
																				   pop.obj_fun.upper_bound(k),
																				   pop.newguys[mutant].genome.parameters(k),
																				   s);

        }
		pop.newguys[mutant].genome.update_chromosomes();
		pop.getFitness(pop.newguys[mutant]);    									//Get Fitness score

		//Perform metropolis
		dH = pop.newguys[mutant].H - pop.guys[mutant].H;		//used to decide if we accept the new guy or not.
		if (dH < 0 || exp(-dH / pop.newguys[mutant].t) > uniform_double_1()) {
			pop.newguys[mutant].born = pop.count.generation;
			pop.copy(pop.guys[mutant], pop.newguys[mutant]);
		}
		else {
			pop.copy(pop.newguys[mutant], pop.guys[mutant]);
		}
	}
}


void crossover(population &pop) {
	int nMatings = (int)(0.2*N);
    Array2d expX;
    Array2d expY;
	int crossoverPoint;
	double rc;
    double dHt0, dHt1;
	double PXX, PYY;			//Selection probabilities
	double PXY = 1, PYX = 1;	//Generating probabilities (PXY = PYX for this operator)
	double TXY, TYX;			//Transition probability
	double ZX, ZY;		        //Sum of Boltzmann-weights for current (X) and offspring(Y) populations

	for (int matings = 0; matings < nMatings; matings++) {
        ZX=0;
        ZY=0;
		pop.roulette_select_two_guys(ZX); //pop.selected 0 will be good, pop.selected 1 random
        expX(0) = exp(-(pop.guys[pop.selected(0)].H-pop.lowest_H)); //good guy
		expX(1) = exp(-(pop.guys[pop.selected(1)].H-pop.lowest_H)); //bad guy
		PXX = ( 1 / ((N - 1)*ZX)*(expX(0) + expX(1))); //P((xi,xj) | x)
		
		//Now mate the newguys to create offsprings
		crossoverPoint = uniform_integer( 1, genomeLength-1);
		for (int i = crossoverPoint; i < genomeLength; i++) {
			pop.newguys[pop.selected(0)].genome.copy_loci(i, pop.guys[pop.selected(1)].genome(i));
			pop.newguys[pop.selected(1)].genome.copy_loci(i, pop.guys[pop.selected(0)].genome(i));
		}

        pop.getFitness(pop.newguys[pop.selected(0)]);
        pop.getFitness(pop.newguys[pop.selected(1)]);

        expY(0) = exp(-(pop.newguys[pop.selected(0)].H-pop.lowest_H));
		expY(1) = exp(-(pop.newguys[pop.selected(1)].H-pop.lowest_H));
		for (int i = 0; i < N; i++) {
			ZY += exp(-(pop.newguys[i].H - pop.lowest_H));
		}
		PYY = (1 / ((N - 1)*ZY)*(expY(0) + expY(1))); //P((yi,yj) | y) selection probability
		dHt0 = (pop.newguys[pop.selected(0)].H - pop.guys[pop.selected(0)].H) / pop.guys[pop.selected(0)].t; //good 
		dHt1 = (pop.newguys[pop.selected(1)].H - pop.guys[pop.selected(1)].H) / pop.guys[pop.selected(1)].t; //bad 

		TXY = PXX*PXY;
		TYX = PYY*PYX;

		rc = exp(-dHt0 - dHt1)*TXY / TYX;
        //Accept or reject
		if (uniform_double_1() < fmin(1, rc)) {
			pop.newguys[pop.selected(0)].born = pop.count.generation;
			pop.newguys[pop.selected(1)].born = pop.count.generation;
			pop.copy(pop.guys[pop.selected(0)], pop.newguys[pop.selected(0)]);
			pop.copy(pop.guys[pop.selected(1)], pop.newguys[pop.selected(1)]);
		}
		else {
			//Revert changes on newguys
			pop.copy(pop.newguys[pop.selected(0)], pop.guys[pop.selected(0)]);
			pop.copy(pop.newguys[pop.selected(1)], pop.guys[pop.selected(1)]);
		}
	}
}



void crossover_elite(population &pop) {
	//Start roulette selection
	int matings;
	int nMatings = (int)(0.2*N);
	int crossoverPoint;
    double ZX,ZY;
	double rc, dHt0, dHt1;
	double PXX, PYY;			//Selection probabilities
	double PXY, PYX;	        //Generating probabilities (PXY = PYX for this operator)
	double TXY, TYX;			//Transition probability
    Array2d expX;
    Array2d expY;
	int random_bestguy1, random_bestguy2; //Guy to inject into pop.selected[1]
	for (matings = 0; matings < nMatings; matings++) {
        ZX = 0;
        ZY = 0;
        //pop.selected 0 is a boltzmann "just good" guy. pop.selected 1 is a random any guy.
        pop.roulette_select_two_guys(ZX); //pop.selected 0 will be good, pop.selected 1 random
        random_bestguy1 = uniform_integer( 0, N_best-1);
		random_bestguy2 = random_bestguy1;
		while(random_bestguy1==random_bestguy2){
			random_bestguy2 = uniform_integer( 0, N_best-1);
		}

        //Make sure the guy pop.selected(0) and random_bestguys aren't already the same guy, or else they will have identical
        //Offspring. Basically we would be copying back a bestguy?
        if ( pop.guys[pop.selected(0)].H == pop.bestguys[random_bestguy1].H ||  pop.guys[pop.selected(0)].H == pop.bestguys[random_bestguy2].H){
            continue;
        }
		if ( pop.guys[pop.selected(1)].H == pop.bestguys[random_bestguy1].H ||  pop.guys[pop.selected(1)].H == pop.bestguys[random_bestguy2].H){
			continue;
		}


        expX(0) = exp(-(pop.guys[pop.selected(0)].H - pop.lowest_H)); //good guy
		expX(1) = exp(-(pop.guys[pop.selected(1)].H - pop.lowest_H)); //gooder guy
        PXX = ( 1 / ((N - 1)*ZX)*(expX(0) + expX(1))); //P((xi,xj) | x)


		//Now mate to create offspring. Let a bestGuy inject DNA in this process!
        //The newguys will be half random bestguy, half guy selected(0)

		bool notequal = pop.newguys[pop.selected(0)].H != pop.newguys[pop.selected(1)].H;
		int num_points = uniform_integer(1, genomeLength/2);
		ArrayXi crossoverPoints(num_points);
		rndChoice(crossoverPoints.data(), num_points,genomeLength-1);
		int mode = EMC_rnd::uniform_integer_1();
		int j = 0;
		for (int i = 0; i < genomeLength;i++){
			if (mode == 0 ){
				pop.newguys[pop.selected(0)].genome.copy_loci(i, pop.guys[pop.selected(1)].genome(i));
				pop.newguys[pop.selected(1)].genome.copy_loci(i, pop.bestguys[random_bestguy1].genome(i));
				if (i == crossoverPoints(j)){ mode = 1; j++;}
			}else if (mode == 1){
				pop.newguys[pop.selected(0)].genome.copy_loci(i, pop.bestguys[random_bestguy1].genome(i));
				pop.newguys[pop.selected(1)].genome.copy_loci(i, pop.guys[pop.selected(0)].genome(i));
				if (i == crossoverPoints(j)){ mode = 0; j++;}
			}

		}

        //From now on the newguys are offsprings. The parents are a good guy pop.selected 1 and a bestguy random_bestguy. Adoptive father is guy pop.selected 1
        pop.getFitness(pop.newguys[pop.selected(0)]);
        pop.getFitness(pop.newguys[pop.selected(1)]);
		expY(0) = exp(-(pop.newguys[pop.selected(0)].H - pop.lowest_H));
		expY(1) = exp(-(pop.newguys[pop.selected(1)].H - pop.lowest_H));
		for (int i = 0; i < N; i++) {
			ZY += exp(-(pop.newguys[i].H - pop.lowest_H));
		}
		PYY = (1 / ((N - 1)*ZY)*(expY(0) + expY(1))); //P((xi,xj) | x)
        PXY = 1;
        PYX = PXY;
		TXY = PXX*PXY;
		TYX = PYY*PYX;

		dHt0 = (pop.newguys[pop.selected(0)].H - pop.guys[pop.selected(0)].H) / pop.guys[pop.selected(0)].t; //good 
		dHt1 = (pop.newguys[pop.selected(1)].H - pop.guys[pop.selected(1)].H) / pop.guys[pop.selected(1)].t; //bad 


		rc = exp(-dHt0 - dHt1)*TXY / TYX;

		//Accept or reject
		if (uniform_double_1() < fmin(1, rc)) {
			//Keep children.
			pop.newguys[pop.selected(0)].born = pop.count.generation;
			pop.newguys[pop.selected(1)].born = pop.count.generation;
			pop.copy(pop.guys[pop.selected[0]], pop.newguys[pop.selected[0]]);
			pop.copy(pop.guys[pop.selected[1]], pop.newguys[pop.selected[1]]);
		}
		else {
			//Keep parents.
			pop.copy(pop.newguys[pop.selected[0]], pop.guys[pop.selected[0]]);
			pop.copy(pop.newguys[pop.selected[1]], pop.guys[pop.selected[1]]);
		}

	}

}

void crossover_smartCopy(population &pop) {
	//Start roulette selection																							
	int matings;
	int nMatings = (int)(0.2*N);
    double ZX, ZY;
    double rc, dHt0, dHt1;
	double PXX, PYY; //Selection probabilities
	double PXY, PYX; //Generating probabilities (PXY != PYX for this operator)
	double TXY, TYX; //Transition probability
    Array2d expX;
    Array2d expY;
	Array4i n_ab_XY; //Exponents for smartCopy probabilities. Note that n_ab.sum() = genomeLength
	Array4i n_ab_YX; //Exponents for smartCopy probabilities. Note that n_ab.sum() = genomeLength
	for (matings = 0; matings < nMatings; matings++) {
        ZX=0;
        ZY=0;
        n_ab_XY.fill(0);
		n_ab_YX.fill(0);
        pop.roulette_select_two_guys(ZX); //pop.selected 0 will be good, pop.selected 1 random
        expX(0) = exp(-(pop.guys[pop.selected(0)].H-pop.lowest_H)); //good guy
        expX(1) = exp(-(pop.guys[pop.selected(1)].H-pop.lowest_H)); //bad guy
        PXX = ( 1 / ((N - 1)*ZX)*(expX(0) + expX(1))); //P((xi,xj) | x)

        //Now mate the newguys to create offsprings
		bitselector_smartCopy(pop, n_ab_XY,n_ab_YX);
		
		//Now we need generating probabilities PXY and PYX
		PXY = 1;
		PYX = 1;
		for (int i = 0; i < 4; i++) {
			PXY *= pow(p_matrix[i], n_ab_XY(i)/genomeLength); //Convert exponent to frequency
			PYX *= pow(p_matrix[i], n_ab_YX(i)/genomeLength); //Convert exponent to frequency
		}

        pop.getFitness(pop.newguys[pop.selected(0)]);
        pop.getFitness(pop.newguys[pop.selected(1)]);

        expY(0) = exp(-(pop.newguys[pop.selected(0)].H-pop.lowest_H));
        expY(1) = exp(-(pop.newguys[pop.selected(1)].H-pop.lowest_H));
        for (int i = 0; i < N; i++) {
            ZY += exp(-(pop.newguys[i].H - pop.lowest_H));
        }
		PYY = (1 / ((N - 1)*ZY)*(expY(0) + expY(1))); //P((yi,yj) | y) selection probability
	    dHt0 = (pop.newguys[pop.selected(0)].H - pop.guys[pop.selected(0)].H) / pop.guys[pop.selected(0)].t; //good
		dHt1 = (pop.newguys[pop.selected(1)].H - pop.guys[pop.selected(1)].H) / pop.guys[pop.selected(1)].t; //bad 
		
		TXY = PXX*PXY;
		TYX = PYY*PYX;

		rc = exp(-dHt0 - dHt1)*TXY / TYX;
		//Accept or reject
		if (uniform_double_1() < fmin(1, rc)) {
			//Keep Children
			pop.newguys[pop.selected(0)].born = pop.count.generation;
			pop.newguys[pop.selected(1)].born = pop.count.generation;
			pop.copy(pop.guys[pop.selected[0]], pop.newguys[pop.selected[0]]);
			pop.copy(pop.guys[pop.selected[1]], pop.newguys[pop.selected[1]]);

		}
		else {
			//Keep parents
			pop.copy(pop.newguys[pop.selected[0]], pop.guys[pop.selected[0]]);
			pop.copy(pop.newguys[pop.selected[1]], pop.guys[pop.selected[1]]);
		}

	}
}

void crossover_snooker(population &pop) {
	//Perhaps the distribution f(r) can be done with a roulette?
	//f(r) = exp(H(x+re))/integral_-r_min^r_max exp(H(x+re)) typ?
	//Notation from snooker crossover
	//Choose 2 guys from pop
	rndChoice(pop.selected.data(), 2, N - 1);
    if (pop.guys[pop.selected(0)].H ==pop.guys[pop.selected(1)].H){
        return;
    }
	pop.line.Through(pop.guys[pop.selected(0)].genome.parameters, pop.guys[pop.selected(1)].genome.parameters);
	//Distance in terms of r between guys
	//Vary r to find the walls of the parameter domain

    double r_max = pop.line.line_max;
    double r_min = pop.line.line_min;
	ArrayXd r_point(r_num);
    if( std::isnan(r_max)){
//		cout << pop.guys[pop.selected(0)].genome.parameters.transpose() << endl;
//		cout << pop.guys[pop.selected(1)].genome.parameters.transpose() << endl;
//        cout << "r_min = "    << r_min << " r_max = " << r_max << endl;
//        cout << "r is nan!! " << endl;
        return;
    }
    if( r_min * r_max > 0){
//		cout << pop.guys[pop.selected(0)].genome.parameters.transpose() << endl;
//		cout << pop.guys[pop.selected(1)].genome.parameters.transpose() << endl;
//        cout << "r_min = "    << r_min << " r_max = " << r_max << endl;
//        cout << "r have the same sign!!! " << endl;
		return;
    }
    if (r_min == r_max){
//		cout << pop.guys[pop.selected(0)].genome.parameters.transpose() << endl;
//		cout << pop.guys[pop.selected(1)].genome.parameters.transpose() << endl;
//        cout << "r_min = "    << r_min << " r_max = " << r_max << endl;
//        cout << "r are equal!!! " << endl << endl;
		return;
    }

    if (fmin(r_min,r_max) > 0.0 ){
//		cout << pop.guys[pop.selected(0)].genome.parameters.transpose() << endl;
//		cout << pop.guys[pop.selected(1)].genome.parameters.transpose() << endl;
//        cout << "r_min = "    << r_min << " r_max = " << r_max << endl;
//        cout << "r_min too large!!! " << endl << endl;
        return;
    }

    if (fmax(r_min,r_max) < 1.0 ){
//		cout << pop.guys[pop.selected(0)].genome.parameters.transpose() << endl;
//		cout << pop.guys[pop.selected(1)].genome.parameters.transpose() << endl;
//
//		cout << "r_min = "    << r_min << " r_max = " << r_max << endl;
//        cout << "r_max too too small!!! "<< endl << endl;
		return;

    }
    if (r_min == 0 || r_max == 0){
//		cout << pop.guys[pop.selected(0)].genome.parameters.transpose() << endl;
//		cout << pop.guys[pop.selected(1)].genome.parameters.transpose() << endl;
//
//        cout << "r_min = "    << r_min << " r_max = " << r_max << endl;
//
//        cout << "r are zero!!! " << endl;

        return;

    }
//    cout << "r_min = "    << r_min << " r_max = " << r_max << endl;
    for (int i = 0; i < r_num; i++) {
		r_point(i) = gaussian_truncated(	//Create gaussians points centered around the good performer
										fmin(r_min, r_max),
										fmax(r_min, r_max),
										1.0,
                                        1.0);
//										0.5*fmax(r_min,r_max));
		pop.snookerGuys[i].genome.set_parameters(pop.line.pointAt(r_point(i)));
		pop.getFitness(pop.line.pointAt(r_point(i)), pop.snookerGuys[i]);
		pop.snookerGuys[i].t = pop.guys[pop.selected(0)].t;
	}


    //Accept or reject
	double lowest_H = pop.snookerGuys[0].H;
	int    lowest_i = 0;
    for (unsigned int i = 0; i < pop.snookerGuys.size(); i++){
        if (pop.snookerGuys[i].H < lowest_H){
			lowest_H = pop.snookerGuys[i].H;
			lowest_i = i;
        }
    }
	double dH = pop.snookerGuys[lowest_i].H - pop.guys[pop.selected(1)].H;
	if (uniform_double_1() < exp(-dH/pop.guys[pop.selected(1)].t)){
		pop.guys[pop.selected(1)].genome.set_parameters(pop.snookerGuys[lowest_i].genome.parameters); //Copy his shit
		pop.guys[pop.selected(1)].H 	= pop.snookerGuys[lowest_i].H;
		pop.guys[pop.selected(1)].born 	= pop.count.generation;
		pop.newguys[pop.selected(1)].born =  pop.count.generation;
        pop.copy(pop.newguys[pop.selected(1)], pop.guys[pop.selected(1)]);
	}





    //Time to make a roulette to see which "r" is chosen.

//	//Make a roulette wheel
//	ArrayXd roulette(r_num);
//	double lucky_number = uniform_double_1();
//
//	double Z = 0;
//    for (int i = 0; i < r_num; i++) {
//		Z += exp(-(pop.snookerGuys[i].H - lowest_H)); //Area on roulette wheel proportional to Boltzmann weights
//		roulette(i) = Z; //Cumulative sum - Low fitness gives large area on roulette
//	}
//	roulette /= Z; //Normalize to a proper distribution
//	for (int i = 0; i < r_num; i++) {
//		if (lucky_number <= roulette(i)) {
//			pop.guys[pop.selected(0)].genome.set_parameters(pop.snookerGuys[i].genome.parameters); //Copy his shit
//			pop.guys[pop.selected(0)].H 	= pop.snookerGuys[i].H;
//			pop.guys[pop.selected(0)].born 	= pop.count.generation;
//			break;
//		}
//	}
}


