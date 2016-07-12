#include<iomanip>
#include<iostream>
#include<fstream>
#include<cmath>
#include<random>
#include<chrono>
#include "memory.h"

#define DELAY_ 1000;
std::ofstream ofile;

double sp_energy(double*&,int);
double tb_energy(double*&,int,double,double);
double total_energy(double*&,int,double,double);

int main(int argc, char* argv[]){
	double 		beta, Z, range = 0.05, V_0 = 1.0,
				ENERGY, energy, new_energy, V1, V2, V1_new, V2_new, delta_energy, prob,
				sig_e=0.0, sig_V1=0.0, sig_V2=0.0, sig_accm= 0.0,
				step_size, pi = acos(-1),
				rand_accept;
	unsigned 	particles, grid, steps, samples, blocks, tb,
				seed_accept, seed_particle, seed_init, seed_move, seed_clock;
	int 		i, j, k, l,
				rand_part, rand_move,
				accepted_moves = 0,
				DELAY = DELAY_;
	char* outfilename;

	if(argc<9){
		std::cout << "Bad usage. Enter also 'particle_count beta tb_0_1 grid steps samples blocks manual_seed_0_1' on same line." << std::endl;
		exit(1);
	}
	else{
		particles = atoi(argv[1]);
		beta = atof(argv[2]);
		tb = atoi(argv[3]);
		grid = atoi(argv[4]);
			step_size = 1.0/(double) grid;
		steps = atoi(argv[5]);
		samples = atoi(argv[6]);
		blocks = atoi(argv[7]);
		seed_clock = atoi(argv[8]);
			if(seed_clock!=0){
				std::cout << "Input initial position seed: ";
				std::cin >> seed_init;
				std::cout << "Input particle seed: ";
				std::cin >> seed_particle;
				std::cout << "Input move seed: ";
				std::cin >> seed_move;
				std::cout << "Input acceptance seed: ";
				std::cin >> seed_accept;
			}
			else{
				seed_init = std::chrono::system_clock::now().time_since_epoch().count();
				for(i=0;i<DELAY;i++){i+1;}
				seed_particle = std::chrono::system_clock::now().time_since_epoch().count();
				for(i=0;i<DELAY;i++){i+1;}
				seed_move = std::chrono::system_clock::now().time_since_epoch().count();
				for(i=0;i<DELAY;i++){i+1;}
				seed_accept = std::chrono::system_clock::now().time_since_epoch().count();
			}
	}

	std::default_random_engine init_generator (seed_init);
	std::default_random_engine part_generator (seed_particle);
	std::default_random_engine move_generator (seed_move);
	std::default_random_engine accept_generator (seed_accept);
	std::uniform_int_distribution<int> position_dist(0,grid-1);
	std::uniform_int_distribution<int> part_dist(0,particles-1);
	std::uniform_int_distribution<int> move_dist(0,1);
	std::uniform_real_distribution<double> accept_dist(0.0,1.0);

	array<double> 	position(particles,0.0),
					new_position(particles,0.0),
					E(blocks,0.0),
					V1_(blocks,0.0),
					V2_(blocks,0.0);
	array<int>		pos(particles,0),
					new_pos(particles,0),
					acc_moves(blocks,0);
	matrix<double>	record(blocks*samples,particles+3,0.0);

	//initialize with random positions
	for(i=0;i<particles;i++){
		pos.memory[i] = position_dist(init_generator);
		new_pos.memory[i] = pos.memory[i];
		position.memory[i] = step_size*(double)pos.memory[i];
		new_position.memory[i] = position.memory[i];
	}

	//initialize energy
	V1 = sp_energy(position.memory,particles);
	V2 = 0.0;
	if(tb!=0){
		V2 = tb_energy(position.memory,particles,range,V_0);
	}
	new_energy = energy = V1+V2;
	ENERGY = 0.0;

	//run algorithm
	//calculate average values after each block
	for(i=0;i<blocks;i++){
		//conduct a set number of samples per block 
		//	evolve if conditions are met
		ENERGY = 0.0;
		for(j=0;j<samples;j++){
			//run for a number of steps before each sample
			for(k=0;k<steps;k++){
				//generator particle and direction
				rand_part = part_dist(part_generator);
				rand_move = position_dist(init_generator);
				//calculate new position
				new_pos.memory[rand_part] = rand_move;
				//calculate new positions given new position indices
				for(l=0;l<particles;l++){
					new_position.memory[l] = step_size*(double)new_pos.memory[l];
				}
				//calculate new energy given new positions
				V1_new = sp_energy(new_position.memory,particles);
				V2_new = 0.0;
				if(tb==1){
					V2_new = tb_energy(new_position.memory,particles,range,V_0);
				}
				new_energy = V1_new + V2_new;
				//calculate delta_e
				delta_energy = new_energy - energy;
				//if delta_e > 0, run metropolis check
				if(delta_energy > 0.0){
					rand_accept = accept_dist(accept_generator);
					prob = exp(-beta*delta_energy);
					//if test is failed, return to old conditions
					if(rand_accept > prob){
						for(l=0;l<particles;l++){
							new_position.memory[l] = position.memory[l];
							new_pos.memory[l] = pos.memory[l];
						}
					}
					//if test is passed, keep new positions, energy, and note accepted move
					else{
						for(l=0;l<particles;l++){
							position.memory[l] = new_position.memory[l];
							pos.memory[l] = new_pos.memory[l];
						}
						V1 = V1_new;
						V2 = V2_new;		
						energy = new_energy;
						accepted_moves++;
					}
					
				}
				//if delta_e <= 0, accept new positions, energy, and note accepted move
				else{
					for(l=0;l<particles;l++){
						position.memory[l] = new_position.memory[l];
						pos.memory[l] = new_pos.memory[l];
					}
					V1 = V1_new;
					V2 = V2_new;
					energy = new_energy;
					accepted_moves++;
				}
			}
			ENERGY += energy;
			V1_.memory[i] += V1;
			V2_.memory[i] += V2;
		}
		E.memory[i] = ENERGY/(double)samples;
		V1_.memory[i] /= (double)samples;
		V2_.memory[i] /= (double)samples;
		acc_moves.memory[i] = accepted_moves;
		accepted_moves = 0;
	}
	ENERGY = 0.0; V1 = 0.0; V2 = 0.0; accepted_moves = 0;
	for(i=0;i<blocks;i++){
		ENERGY += E.memory[i];
		V1 += V1_.memory[i];
		V2 += V2_.memory[i];
		accepted_moves += acc_moves.memory[i];
	}
	ENERGY /= (double)blocks;
	V1 /= (double)blocks;
	V2 /= (double)blocks;

	for(i=0;i<blocks;i++){
		sig_e = fabs((E.memory[i] - ENERGY)*(E.memory[i] - ENERGY));
		sig_V1 = fabs((V1_.memory[i] - V2)*(V1_.memory[i] - V1));
		sig_V2 = fabs((V2_.memory[i] - V2)*(V2_.memory[i] - V2));
		sig_accm = fabs((acc_moves.memory[i] - accepted_moves)*(acc_moves.memory[i] - accepted_moves));
	}
	sig_e = pow(sig_e,0.5)/(double) blocks;
	sig_V1 = pow(sig_V1,0.5)/(double) blocks;
	sig_V2 = pow(sig_V2,0.5)/(double) blocks;
	sig_accm = pow(sig_accm,0.5)/(double) (blocks*blocks);

	std::cout << "energy: " << ENERGY << " (" << sig_e << ")" << std::endl;
	std::cout << "<V1>: " << V1 << " (" << sig_V1 << ")" << std::endl;
	std::cout << "<V2>: " << V2 <<  " (" << sig_V2 << ")" <<std::endl;
	std::cout << "<AM>: " << (double) accepted_moves/(double)blocks << " (" << sig_accm << ")" << std::endl;
				
	return 0;
}

double sp_energy(double*& position, int particles){
	int 	i;
	double 	energy = 0.0,
			pi = acos(-1);
	for(i=0;i<particles;i++){
		energy += cos(20*pi*position[i]);
	}
	energy += (double)particles;
	return energy;
}
double tb_energy(double*& position, int particles, double range, double V_0){
	int 	i, j;
	double 	energy = 0.0,
			distance;
	for(i=0;i<particles;i++){
		for(j=i+1;j<particles;j++){
			distance = fabs(position[i] - position[j]);
			if(distance < range){
				energy += V_0;
			}
		}
	}
	return energy;
}
double total_energy(double*& position, int particles, double range, double V_0){
	int 	i, j;
	double 	energy = 0.0,
			distance,
			pi = acos(-1);
	for(i=0;i<particles;i++){
		for(j=i+1;j<particles;j++){
			distance = fabs(position[i] - position[j]);
			if(distance < range){
				energy += V_0;
			}
		}
		energy += cos(20*pi*position[i]);
	}
	energy += (double) particles;
	return energy;
}
