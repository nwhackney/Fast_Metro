#include "site.hpp"
#include <vector>
#include <math.h>
#include <iostream>

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

	void init(int Lattice_size, int Number);
	void set_const(double j, double sig, double frustration);
	void restart(int Length, int Number, std::string infile);

	void print_data(std::string file_name);
	void print_gnu(std::string file_name);

	void Metropolis(double T, std::ofstream &Efile, std::vector<double> &accepted,double width);
	void check();

	void rotate(int i, double theta);
	void flip(int i);

	int occupied(int i);
	int how_many();

	double angle(int i);
	double H_local(int i);
	double H();
	double strain(int i);

};