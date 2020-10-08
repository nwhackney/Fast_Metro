#include <iostream>
#include <fstream>
#include <sstream>

#include "../include/lattice.hpp"
#include "../include/simulation.hpp"
#include "../include/tinytoml-master/include/toml/toml.h"
#include "../include/gsl/gsl_rng.h"

using namespace std;

void simulation::read_config()
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
	const toml::Value* Sp = v.find("Sigma");
	const toml::Value* fp = v.find("f");
	const toml::Value* Tp = v.find("Time");
	const toml::Value* Sl = v.find("Slope");
	const toml::Value* ET = v.find("End_Time");
	const toml::Value* Etemp = v.find("End_Temp");
	
	const toml::Value* outp = v.find("output");
	
	const toml::Value*  W = v.find("Width");
	
	const toml::Value* rst = v.find("Restart_Time");
	const toml::Value* im = v.find("Restart");
	const toml::Value* i = v.find("In_File");
	const toml::Value* s=v.find("Seed");
	const toml::Value* nr=v.find("Number_of_Replicas");
	const toml::Value* sp=v.find("Swap_Probability");
	const toml::Value* rt=v.find("Run_Type");
	
	
	L= Np->as<int>();
	N= Occp->as<int>();
	Time=Tp->as<int>();
	end_time=ET->as<int>();
	restart_t=rst->as<int>();
	seed=s->as<int>();
	Num_Rep=nr->as<int>();
	
	slp=Sl->as<double>();
	w=W->as<double>();
	J=Jp->as<double>();
	double Sigma=Sp->as<double>();
	K=Sigma-J;
	f=fp->as<double>();
	Tf=Etemp->as<double>();
	SP=sp->as<double>();

	restart=im->as<string>();
	out_file=outp->as<string>();
	in_file=i->as<string>();
	run_type=rt->as<string>();
}

void simulation::simulated_annealing()
{
	//initializing gsl random number generator
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
		Temp=(1.0/cosh(w*slope*((double) t)))+Tf;
		
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
	inf<<"Gamma= "<<(J+K)/abs(J)<<endl;
	inf<<"Final Energy: "<<crystal.H()<<endl;
	inf<<"Time: "<<duration<<endl;

	inf.close();
	Edat.close();
}

void simulation::parallel_tempering()
{
	const gsl_rng_type *T;
	gsl_rng_env_setup();
	T=gsl_rng_default; // Currently this is the default Tausworthe engine, Chris suggested the Mersenee Twister. Should maybe look into switching/why
	r=gsl_rng_alloc(T);

	gsl_rng_set(r,seed);

	vector<int> seeds;
	seeds.resize(Num_Rep);
	for (int i=0; i<Num_Rep; i++)
	{seeds[i]=gsl_rng_get(r);}

	lattice crystal;
	crystal.set_const(J,K,f);

	crystal.init(L,N,r);

	crystal.print_data("init_data.dat");
	crystal.print_gnu("init");

	// Initializing Ensemble

	vector<lattice> ensemble;
	ensemble.resize(Num_Rep,crystal);

	vector<double> temp;
	temp.resize(Num_Rep);

	temp[0]=0.0;
	for (int i=0; i<Num_Rep; i++)
	{
		temp[i]= 1.0/(Num_Rep-i); // Geometric Temperature Distribution. Should Change to Iterative method. See Rathore et. al. for details
	}

	vector<gsl_rng *> rng;
	rng.resize(Num_Rep);

	gsl_rng_env_setup();
	T=gsl_rng_default; // Currently this is the default Tausworthe engine, Chris suggested the Mersenee Twister. Should maybe look into switching/why
	
	int inc=0;
	for (auto &i: rng)
	{
		i=gsl_rng_alloc(T);
		gsl_rng_set(i,seeds[inc]);
		inc++;
	}

	cout<<"Initial Energy: "<<crystal.H()<<endl;

	stringstream out;
	out<<out_file;

	vector<ofstream> Efiles;
	Efiles.resize(Num_Rep);
	for (int i=0; i<Num_Rep; i++)
	{
		stringstream nom;
		nom<<"Energy_"<<temp[i]<<".dat";
		Efiles[i].open(nom.str());
	}

	double duration;
	clock_t start;
		start = clock();

	vector <vector<double>> accepted;
	accepted.resize(Num_Rep);
	for (auto &i: accepted) {i.resize(8,0.0);}

	ofstream rot_acc, tran_acc, glide_acc, loc_acc;
	rot_acc.open("Rotation_Acceptance.dat");
	tran_acc.open("Translation_Acceptance.dat");
	glide_acc.open("Glide_Acceptance.dat");
	loc_acc.open("Local_Acceptance.dat");

	// Parallel Tempering Starts Here

	double swap_accpeted=0.0,
		  swap_tried=0.0;

	for (int time=0; time<=Time; time++)
	{
		for (int n=0; n<Num_Rep; n++)
		{
			ensemble[n].Metropolis(temp[n],Efiles[n],accepted[n],rng[n]);
		
			if (n==Num_Rep-1) {continue;} //Skips trying to swap highest temperature replica, as there is nothing to swap with...
			double flag=gsl_rng_uniform(r); //Now the replica swaps happens randomly (after lattice sweeps) instead of on a schedule. Suggested by Falcioni & Deem
			if (flag <= SP)
			{
				swap_tried+=1.0;
				double E_i, E_j, b_i, b_j, expo, alpha, U;
				lattice temporary;
				
				E_i=ensemble[n].H();
				b_i=1.0/temp[n];

				E_j=ensemble[n+1].H();
				b_j=1.0/temp[n+1];

				expo=(b_j-b_i)*(E_j-E_i);

				alpha=gsl_rng_uniform(r);

				U= exp(expo);
				if (alpha < fmin(1.0,U))
				{
					swap_accpeted+=1.0;
					temporary=ensemble[n];
					ensemble[n]=ensemble[n+1];
					ensemble[n+1]=temporary;	
				}
			}
		}
	}

	cout<<swap_accpeted/swap_tried<<endl;

	for (int t=0; t<10000; t++) // (Apparently) Running lowest temperature (T=0) ensemble for additional timesteps 
	{
		ensemble[0].Metropolis(temp[0],Efiles[0],accepted[0],rng[0]);
	}

	// Parallel Tempering Ends Here

	gsl_rng_free(r); //Freeing the rng
	for (auto &i: rng){gsl_rng_free(i);}

	for (auto &i: Efiles) {i.close();}

	rot_acc.close(); tran_acc.close(); glide_acc.close(); loc_acc.close();

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	cout<<"Time: "<<duration<<endl;
	cout<<"Final Energy: "<<ensemble[0].H()<<endl;

	ensemble[0].print_data(out.str());
	ensemble[0].print_gnu(out.str());

	stringstream info;
	info<<"info_"<<out_file<<".dat";

	ofstream inf;
	inf.open(info.str());
	inf<<L<<"x"<<L<<" lattice"<<endl;
	inf<<N<<" spin sites"<<endl;
	inf<<Time<<" sweeps"<<endl;
	inf<<"J="<<J<<" K="<<K<<" f="<<f<<endl;
	inf<<"Gamma= "<<(J+K)/abs(J)<<endl;
	inf<<"Final Energy: "<<ensemble[0].H()<<endl;
	inf<<"Time: "<<duration<<endl;

	inf.close();
}

void simulation::run()
{
	if (run_type=="Simulated_Annealing")
	{
		simulated_annealing();
	}
	else if (run_type=="Parallel_Tempering")
	{
		parallel_tempering();
	}

}