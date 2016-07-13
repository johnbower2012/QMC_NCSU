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
	int 	i, j, k, l,
			steps, walkers0, walkers, walkers_t,
			accm=0, mult;
	double	alpha, omega,
			sigma, random, etrial,
			energy, energy_prime,
			x, x_prime,
			weights, tau=0.0, delta_tau, delta_tau_sqrt;
	char* filename;

	if(argc<8){
		std::cout << "Bad usage. Enter also 'omega alpha etrial steps delta_tau walkers filename' on same line." << std::endl;
		exit(1);
	}
	else{
		omega = atof(argv[1]);
		alpha = atof(argv[2]);
		etrial = atof(argv[3]);
		steps = atoi(argv[4]);
		delta_tau = atof(argv[5]);
			delta_tau_sqrt = pow(delta_tau,0.5);
		walkers = walkers0 = walkers_t = atoi(argv[6]);
		filename = argv[7];
	}

	unsigned seed_n = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator_n(seed_n);
	std::uniform_real_distribution<double> real_dist(0.0,1.0);
	std::normal_distribution<double> norm_dist(0.0,1.0);
	unsigned seed_r = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator_r(seed_r);


	matrix<double> 	walker(steps,walkers,0.0),
					walker_t(steps,walkers,0.0);
	array<double>	ENERGY(walkers,0.0),
					energy_record(steps,0.0),
					weight(walkers,0);
	array<int> 		multiplicity(walkers,0);

	for(i=0;i<walkers;i++){
		random = norm_dist(generator_n);
		walker.memory[0][i] = x = random;
		ENERGY.memory[i] = gaussian_energy(alpha,omega,x);
		energy_record.memory[0] += ENERGY.memory[i];
	}
	energy_record.memory[i] /= walkers;

	for(i=1;i<steps;i++){
		//loop over walkers
		for(j=0;j<walkers;j++){
			x = walker.memory[i-1][j];
			//energy = gaussian_energy(alpha,omega,x);

			random = norm_dist(generator_n);
			x_prime = x - delta_tau*alpha*x + random*delta_tau_sqrt;
			walker.memory[i][j] = x_prime;
			energy_prime = gaussian_energy(alpha,omega,x_prime);

			weight.memory[j] = exp(-(energy_prime - etrial)*delta_tau);
		}

		//count new multiplicity of each state
		//weight towards original count
		walkers_t = 0;
		for(j=0;j<walkers;j++){
			weight.memory[j] *= (double) walkers0/(double) walkers;
			multiplicity.memory[j] = weight.memory[j] + real_dist(generator_r);
			walkers_t += multiplicity.memory[j];
		}
		
		//BRANCHING

		//walker count test
		if(walkers_t==0){
			std::cout << "All walkers died on step " << i << "." << std::endl;
			break;
		}
		else if(walkers_t>2*walkers0){
			std::cout << "Walkers exceeded " << 2*walkers0 << " on step " << i << "." << std::endl;
			break;
		}
		//construct new branch set		
		walker_t.resize(steps,walkers_t,false);
		for(j=0;j<i+1;j++){
			l=0;
			for(k=0;k<walkers;k++){
				mult = multiplicity.memory[k];
				while(mult>0){
					walker_t.memory[j][l] = walker.memory[j][k];
					mult -= 1;
					l++;
				}
			}
		}
		walkers = walkers_t;
		walker.resize(steps,walkers,false);
		multiplicity.resize(walkers,false);
		ENERGY.resize(walkers,false);
		for(k=0;k<walkers;k++){
			for(j=0;j<i+1;j++){
				walker.memory[j][k] = walker_t.memory[j][k];
			}
		}

		//energies
		for(j=0;j<walkers;j++){
			x = walker.memory[i][j];
			ENERGY.memory[j] = gaussian_energy(alpha,omega,x);
			energy_record.memory[i] += ENERGY.memory[j];
		}
		energy_record.memory[i] /= walkers;

		//tau
		tau += delta_tau;
	}


	energy = 0.0;
	for(i=0;i<walkers;i++){
		energy += ENERGY.memory[i];
	}
	energy /= walkers;
	ofile.open(filename);	
	for(i=1;i<steps;i*=10){
		ofile << delta_tau*(i+1) << std::setw(15) << energy_record.memory[i];
		for(j=0;j<walkers;j++){
			ofile << std::setw(15) << walker.memory[i][j];
		}
		ofile << std::endl;
	}
	ofile.close();

	std::cout << "walkers: " << walkers << std::endl;
	std::cout << "energy: " << energy << std::endl;
	std::cout << "tau: " << tau << std::endl;

	return 0;
}
