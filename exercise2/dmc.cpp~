#include<iomanip>
#include<iostream>
#include<fstream>
#include<cmath>
#include<chrono>
#include<random>
#include "memory.h"

std::ofstream ofile;

inline double gaussian_psi(double alpha, double x){
	return exp(-alpha*x*x/2.0);
}
inline double gaussian_energy(double alpha, double omega, double x){
	return -0.5*(alpha*alpha*x*x - alpha) + 0.5*omega*omega*x*x;
}

int main(int argc, char* argv[]){
	unsigned seed;
	int 	i, j, k, l,
			blocks, samples, run, sample_walkers,
			walkers0, walkers, walkers_t,
			accm=0, mult;
	double	alpha, omega,
			sigma, random, etrial,
			energy_prime,
			x, x_prime,
			E=0.0, E_sig=0.0,
			weights, tau=0.0, delta_tau, delta_tau_sqrt;
	char* filename;

	if(argc<11){
		std::cout << "Bad usage. Enter also 'omega alpha etrial blocks steps run delta_tau walkers seed filename' on same line." << std::endl;
		exit(1);
	}
	else{
		omega = atof(argv[1]);
		alpha = atof(argv[2]);
		etrial = atof(argv[3]);
		blocks = atoi(argv[4]);
		samples = atoi(argv[5]);
		run = atoi(argv[6]);
		delta_tau = atof(argv[7]);
		walkers = walkers0 = walkers_t = atoi(argv[8]);
		seed = atoi(argv[9]);
		filename = argv[10];

		sample_walkers = samples*walkers;
		delta_tau_sqrt = pow(delta_tau,0.5);
		if(seed==0){
			seed = std::chrono::system_clock::now().time_since_epoch().count();
		}	
	}


	std::uniform_real_distribution<double> real_dist(0.0,1.0);
	std::normal_distribution<double> norm_dist(0.0,1.0);
	std::default_random_engine generator(seed);

	array<double>	walker(walkers,0.0);
	matrix<double>	energy(walkers,2,0.0),
					energy_block(blocks,2,0.0);

	for(i=0;i<walkers;i++){
		walker.memory[i] = norm_dist(generator);
	}

	for(i=0;i<blocks;i++){
		//reset energies
		for(j=0;j<walkers;j++){
			energy.memory[j][0] = 0.0;
			energy.memory[j][1] = 0.0;
		}
		for(j=0;j<samples;j++){
			for(k=0;k<walkers;k++){
				for(l=0;l<run;l++){
					x = walker.memory[k];
					random = norm_dist(generator);
					x_prime = x - delta_tau*alpha*x + random*delta_tau_sqrt;
					energy_prime = gaussian_energy(alpha,omega,x_prime);
					walker.memory[k] = x_prime;
				}

				weights = exp(-(energy_prime - etrial)*delta_tau);
				energy.memory[k][0] = energy_prime*weights;
				energy.memory[k][1] = energy.memory[k][0]*energy.memory[k][0];
			}
			for(k=0;k<walkers;k++){
				energy_block.memory[i][0] += energy.memory[k][0];
				energy_block.memory[i][1] += energy.memory[k][1];
			}
		}
		energy_block.memory[i][0] /= (double) (sample_walkers);
		energy_block.memory[i][1] /= (double) (sample_walkers);
		E += energy_block.memory[i][0];
		E_sig += energy_block.memory[i][1];
	}
	E /= (double) blocks;
	E_sig /= (double) blocks; 
	E_sig = pow(fabs(E_sig - E*E),0.5);

	std::cout << E << std::setw(15) << E_sig << std::endl;


	return 0;
}
