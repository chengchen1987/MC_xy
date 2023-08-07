#pragma once

#include <random>

class Square_Lattice
{
private:
	int Lx, Ly;
public:
	Square_Lattice(int lx, int ly);

	int get_Lx(void) { return Lx; }
	int get_Ly(void) { return Ly; }
	void NN_sites(int site, int* nn_sites);
	void print_spin(double* spins);
};

// we print basic measurements of each 'bin', and post-compute observables   
class MC_Measurements
{
public:
	double E, E2, M, M2, M4;
};

class MC_xy
{
private:
	Square_Lattice* lattice;

public:
	MC_xy(Square_Lattice* _lattice);
	~MC_xy();

	// lattice information 
	int Lx, Ly, N;

	// MC parameters
	int n_warmup, n_measure, n_bins;
	double T_min, T_inc, n_T;
	double* T_vec;
	
	// MC functions 
	void MC_warmup(int iT, double T, double* spin);
	
	void MC_measure(int iT, double T, double* spin);
	
	void MC_step_measure(double T, double *spin, MC_Measurements &obs);
	double MC_step_Metropolis(double T, double* spin);
	double calc_H(double* spin);
	double calc_M(double* spin);
	void cal_mean_measurements(int dlen, MC_Measurements *step_obs, MC_Measurements &obs);

	std::random_device rd;  // Will be used to obtain a seed for the random number engine
};
