#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<random>
#include<chrono>
#include "memory.h"

#define RANGE 5.0

std::ofstream ofile;

inline double gaussian_psi(double alpha, double x){
	return exp(-alpha*x*x/2.0);
}
inline double gaussian_energy(double alpha, double omega, double x){
	return -0.5*(alpha*alpha*x*x - alpha) + 0.5*omega*omega*x*x;
}

int main(int argc, char* argv[]){
	//declare variables
	int 	i, j, k, l,
			steps, walkers, run,
			accm=0;
	double 	alpha, omega,
			sigma, random,
			psi, psi_new, prob,
			x, x_new,
			ENERGY_;
	char* filename;

	if(argc<6){
		std::cout << "Bad usage. Enter also 'steps walkers run sigma filename' on same line." << std::endl;
		exit(1);
	}
	else{
		steps = atoi(argv[1]);
		walkers = atoi(argv[2]);
		run = atoi(argv[3]);
		sigma = atof(argv[4]);
		filename = argv[5];
	}

	//random generation
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::uniform_real_distribution<double> real_dist(0.0,1.0);
	std::default_random_engine generator(seed);

	//declare arrays
	matrix<double> 	walker(steps,walkers,0.0),
					energy(steps,walkers,0.0);
	array<double>	ENERGY(walkers,0.0);

	//initialize system	
	for(i=0;i<walkers;i++){
		random = real_dist(generator);
		walker.memory[0][i] = x = (random - 0.5)*RANGE*2.0;
		energy.memory[0][i] = gaussian_energy(alpha,omega,x);
	}

	ofile.open(filename);
	omega = 1;
	for(l=1;l<20;l++){
		//initialize system	
		alpha = (double) l/5.0;
		accm = 0; ENERGY_ = 0.0;
		for(i=0;i<walkers;i++){
			random = real_dist(generator);
			walker.memory[0][i] = x = (random - 0.5)*RANGE*2.0;
			energy.memory[0][i] = gaussian_energy(alpha,omega,x);
			ENERGY.memory[i] = 0.0;
		}
		//run MC
		for(i=1;i<steps;i++){
			for(j=0;j<walkers;j++){
				x = walker.memory[i-1][j];
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
				energy.memory[i][j] = gaussian_energy(alpha,omega,x);
				walker.memory[i][j] = x;
			}
		}
		for(i=0;i<walkers;i++){
			for(j=0;j<steps;j++){
				ENERGY.memory[i] += energy.memory[j][i];
			}
			ENERGY.memory[i] /= steps;
			ENERGY_ += ENERGY.memory[i];
		}
		ENERGY_ /= walkers;
		ofile << alpha << std::setw(15) << ENERGY_ << std::setw(15) << 0.25*(omega*omega - alpha*alpha)/alpha + alpha/2.0 << std::setw(15) << (double)accm/(double)(run*steps*walkers) << std::endl;
	}
	ofile.close();


	return 0;
}
