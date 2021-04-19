#include "site.hpp"
#include <vector>
#include <math.h>
#include <iostream>
#include "../include/gsl/gsl_rng.h"

class lattice
{

	int L,
	    N,
	    V;

	double J,K,f;

	std::vector<site> lattice;
	std::vector<int> occ,
	                 vac;

public:

	void init(int Lattice_size, int Number, gsl_rng * rng);
	void square_init(int Lattice_size, int W, gsl_rng * rng);
	void rect_init(int Lx, int Ly, int W, gsl_rng * rng);
	void ribbon_init(int L, int width, gsl_rng * rng);
	void set_const(double j, double sig, double frustration);
	void restart(int Length, int Number, std::string infile);

	void print_data(std::string file_name);
	void print_gnu(std::string file_name);

	void spin_correlation(std::ofstream &scfile);

	void Metropolis(double T, std::ofstream &Efile, std::vector<double> &accepted, gsl_rng * rng);
	void Spin_Metropolis(double T, std::ofstream &Efile, std::vector<double> &accepted, gsl_rng * rng);
	void Local_Metropolis(double T, std::ofstream &Efile, std::vector<double> &accepted, gsl_rng * rng);
	void check();

	void rotate(int i, double theta);
	void flip(int i);

	int occupied(int i);
	int how_many();
	int size();
	int occ_to_index(int i);

	double angle(int i);
	double H_local(int i);
	double H();
	double strain(int i);

};