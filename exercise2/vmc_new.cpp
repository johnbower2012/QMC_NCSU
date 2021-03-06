#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<random>
#include<chrono>
#include "memory.h"

#define RANGE 5.0

std::ofstream ofile;

inline double gaussian_psi(double*& param){
	return exp(-param[1]*param[0]*param[0]/2.0);
}
inline double gaussian_energy(double*& param){
	return -0.5*(param[1]*param[1]*param[0]*param[0] - param[1]) + 0.5*param[2]*param[2]*param[1]*param[1];
}	
/*************************************
	coord[0] = CURRENT position
	coord[1] = NEXT position
	coord[2] = random step width

	param used as input for function
	param[0] must be coordinate
*************************************/
template <typename FUNC>
bool metropolis(double*& coord, std::default_random_engine generator,  double*& param, FUNC function){
	double 	psi1, psi2, x2, ratio,
			x1 = coord[0],
			sigma = coord[2];

	std::uniform_real_distribution real_dist(0.0,1.0);

	//current values
	psi1 = function(param);

	//new values
	random = real_dist(generator);
	x2 = x1 + sigma*(random - 0.5);
	param[0] = x2;
	psi2 = function(param);

	//comparison
	if(psi2 < psi1){
		ratio = psi2/psi1;
		random = real_dist(generator);
		//deny
		if(random > ratio){
		}
		//accept
		else{
			x1 = x2;
		}
	}
	//accept
	else{
		x1= x2;
	}
}


int main(int argc, char* argv[]){
	//declare variables
	int 	i, j, k, l, m,
			blocks, samples, walkers, run,
			accm=0, totm=0, sample_walkers=0;
	double 	alpha, omega,
			sigma, random,
			psi, psi_new, prob,
			x, x_new,
			E, E_SQ, accmratio;
	char* filename;

	if(argc<7){
		std::cout << "Bad usage. Enter also 'blocks samples run walkers sigma filename' on same line." << std::endl;
		exit(1);
	}
	else{
		blocks = atoi(argv[1]);
		samples = atoi(argv[2]);
		run = atoi(argv[3]);
		walkers = atoi(argv[4]);
		sigma = atof(argv[5]);
		filename = argv[6];
		omega = 1;
		totm = samples*run*walkers;
		sample_walkers = samples*walkers;
	}

	//random generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::uniform_real_distribution<double> real_dist(0.0,1.0);
	std::default_random_engine generator(seed);

	//declare arrays
	array<double>	energy(walkers,0.0),
					energy_sq(walkers,0.0),
					walker(walkers,0.0),
					acc_move(walkers,0.0),
					acc_move_block(blocks,0.0);
	matrix<double>	energy_block(blocks,2,0.0);

	//initialize system
	for(i=0;i<walkers;i++){
		random = real_dist(generator);
		walker.memory[i] = x = (random - 0.5)*RANGE*2.0;
		energy.memory[i] = gaussian_energy(alpha,omega,x);
	}
	//open file
	ofile.open(filename);
	for(i=0;i<walkers;i++){
		random = real_dist(generator);
		walker.memory[i] = (random - 0.5)*RANGE*2.0;
	}
	for(m=1;m<20;m++){
		//initialize system	
		alpha = (double) m/5.0;
		//run MC
		for(l=0;l<blocks;l++){
			//reset energies
			accm = 0;
			for(i=0;i<walkers;i++){
				energy.memory[i] = 0.0;
				energy_sq.memory[i] = 0.0;
			}
			//samples within each block
			for(i=0;i<samples;i++){
				//loop over walkers
				for(j=0;j<walkers;j++){
					//set initial x coordinate
					x = walker.memory[j];
					//loop over number of runs before statistics are calculated
					for(k=0;k<run;k++){
						//current values
						psi = pow(gaussian_psi(alpha,x),2);

						//new values
						random = real_dist(generator);
						x_new = x + sigma*(random - 0.5);
						psi_new = pow(gaussian_psi(alpha,x_new),2);

						//comparison
						if(psi_new < psi){
							prob = psi_new/psi;
							random = real_dist(generator);
							//deny
							if(random > prob){
							}
							//accept
							else{
								x = x_new;
								accm++;
							}
						}
						//accept
						else{
							x = x_new;
							accm++;
						}
					}
					//update walker position && energy statistics
					walker.memory[j] = x;
					energy.memory[j] += gaussian_energy(alpha,omega,x);
					energy_sq.memory[j] += energy.memory[j]*energy.memory[j];
				}
			}
			//block statistics
			for(i=0;i<walkers;i++){
				energy_block.memory[l][0] += energy.memory[i];
				energy_block.memory[l][1] += energy_sq.memory[i];
			}
			energy_block.memory[l][0] /= (double)sample_walkers;
			energy_block.memory[l][1] /= (double)(sample_walkers*sample_walkers);
			acc_move_block.memory[l] = (double)accm/(double)totm;

			E = energy_block.memory[l][0];
			E_SQ = energy_block.memory[l][1];
			energy_block.memory[l][1] = pow(E_SQ - E*E,0.5);
		}
		//average block statistics
		E = E_SQ = accmratio = 0.0;
		for(l=0;l<blocks;l++){
			E += energy_block.memory[l][0];
			E_SQ += energy_block.memory[l][1];
			accmratio += acc_move_block.memory[l];
		}
		E /= (double)blocks;
		E_SQ /= (double)blocks;
		accmratio /= (double)blocks;
		//write to filename
		ofile << alpha << std::setw(15) << 0.25*(omega*omega - alpha*alpha)/alpha + alpha/2.0 << std::setw(15) << E << std::setw(15) << E_SQ << std::setw(15) << accmratio << std::endl;
	}
	ofile.close();


	return 0;
}
