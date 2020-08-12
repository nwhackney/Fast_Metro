#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <limits>

#include "gsl/gsl_rng.h"
#include "../include/lattice.hpp"
#include "../include/tinytoml-master/include/toml/toml.h"

using namespace std;

void run_config()
{
	std::ifstream ifs("config.toml");
	toml::ParseResult pr = toml::parse(ifs);

	if (!pr.valid())
	{
    		cout << pr.errorReason << endl;
     	return;
	}

	const toml::Value& v = pr.value;
	const toml::Value* Np = v.find("L");
	const toml::Value* Occp = v.find("N");
	const toml::Value* Jp = v.find("J");
	const toml::Value* Kp = v.find("K");
	const toml::Value* fp = v.find("f");
	const toml::Value* Tp = v.find("Time");
	const toml::Value* Sl = v.find("Slope");
	const toml::Value* ET = v.find("End_Time");
	
	const toml::Value* outp = v.find("output");
	
	const toml::Value*  W = v.find("Width");
	
	const toml::Value* rst = v.find("Restart_Time");
	const toml::Value* im = v.find("Restart");
	const toml::Value* i = v.find("In_File");
	const toml::Value* s=v.find("Seed");
	
	
	int L= Np->as<int>();
	int N= Occp->as<int>();
	int Time=Tp->as<int>();
	int end_time=ET->as<int>();
	int restart_t=rst->as<int>();
	int seed=s->as<int>();
	
	double slp=Sl->as<double>();
	double w=W->as<double>();
	double J=Jp->as<double>();
	double K=Kp->as<double>();
	double f=fp->as<double>();
	string restart=im->as<string>();
	string out_file=outp->as<string>();
	string in_file=i->as<string>();

	//initializing gsl random number generator
	gsl_rng * r;
	const gsl_rng_type *T;
	gsl_rng_env_setup();
	T=gsl_rng_default; // Currently this is the default Tausworthe engine, Chris suggested the Mersenee Twister. Should maybe look into switching/why
	r=gsl_rng_alloc(T);

	gsl_rng_set(r,seed);

	lattice crystal;
	crystal.set_const(J,K,f);

	if (restart=="yes")
	{
		crystal.restart(L,N,in_file);
		//crystal.print_data("init");
	}
	else
	{
		crystal.init(L,N,r);
		crystal.print_data("init_data.dat");	
	}

	crystal.print_gnu("init");

	cout<<"Initial Energy: "<<crystal.H()<<endl;

	stringstream Efile;
	Efile<<"Energy_"<<out_file<<".dat";
	stringstream out;
	out<<out_file;

	ofstream Edat;
	Edat.open(Efile.str());

	double duration;
	clock_t start;
		start = clock();
	
	double slope,
	       Temp;

	vector<double> accepted;
	accepted.resize(8,0.0);

	ofstream rot_acc, tran_acc, glide_acc, loc_acc;
	rot_acc.open("Rotation_Acceptance.dat");
	tran_acc.open("Translation_Acceptance.dat");
	glide_acc.open("Glide_Acceptance.dat");
	loc_acc.open("Local_Acceptance.dat");

	for (int t=restart_t; t<Time; t++)
	{
		slope=10.0/(slp);
		Temp=1.0/cosh(w*slope*((double) t));
		
		crystal.Metropolis(Temp,Edat,accepted, r);

		// if (t%100000==0)
		// {
		// 	stringstream look;
		// 	look<<"snap_"<<t;
		// 	crystal.print_gnu(look.str());
		// }

		if (t%1000==0)
		{
			rot_acc<<t<<" "<<accepted[1]/accepted[0]<<endl;
			tran_acc<<t<<" "<<accepted[3]/accepted[2]<<endl;
			glide_acc<<t<<" "<<accepted[5]/accepted[4]<<endl;
			loc_acc<<t<<" "<<accepted[7]/accepted[6]<<endl;
		}

		if (t%50000==0)
		{
			stringstream inter;
			inter<<"Sys";
			crystal.print_data(inter.str());
		}
	}

	for (int t=0; t<=end_time; t++)
	{
		crystal.Metropolis(0.0,Edat,accepted, r);
	}

	gsl_rng_free(r); //Freeing the rng

	rot_acc.close(); tran_acc.close(); glide_acc.close(); loc_acc.close();

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	cout<<"Time: "<<duration<<endl;
	cout<<"Final Energy: "<<crystal.H()<<endl;

	crystal.print_data(out.str());
	crystal.print_gnu(out.str());

	stringstream info;
	info<<"info_"<<out_file<<".dat";

	ofstream inf;
	inf.open(info.str());
	inf<<L<<"x"<<L<<" lattice"<<endl;
	inf<<N<<" spin sites"<<endl;
	inf<<Time<<" sweeps"<<endl;
	inf<<"J="<<J<<" K="<<K<<" f="<<f<<endl;
	inf<<"Gamma= "<<(abs(J)-K)/abs(J)<<endl;
	inf<<"Final Energy: "<<crystal.H()<<endl;
	inf<<"Time: "<<duration<<endl;

	inf.close();
	Edat.close();

}

int main()
{
	run_config();

	return 0;
}