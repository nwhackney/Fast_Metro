#include <iostream>
#include <fstream>
#include <sstream>

#include "../include/lattice.hpp"
#include "../include/simulation.hpp"
#include "../include/Hoshen_Kopelman.hpp"
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
	Sigma=Sp->as<double>();
	K=Sigma-J;
	f=fp->as<double>();
	Tf=Etemp->as<double>();
	SP=sp->as<double>();

	restart=im->as<string>();
	out_file=outp->as<string>();
	in_file=i->as<string>();
	run_type=rt->as<string>();
}

float clamp(float x, float lowerlimit, float upperlimit)
{
  if (x < lowerlimit)
    x = lowerlimit;
  if (x > upperlimit)
    x = upperlimit;
  return x;
}

double Step_Temp(double t)
{
	double T=0.0;
	if (0.0<=t and t<150000.0){T=1.0;}
	else if (200000.0<=t and t<350000.0){T=0.8;}
	else if (400000.0<=t and t<550000.0){T=0.6;}
	else if (600000.0<=t and t<750000.0){T=0.4;}
	else if (800000.0<=t and t<950000.0){T=0.2;}
	else if (1000000.0<=t and t<=1500000.0){T=0.0;}

	else if (150000.0<=t and t<200000.0)
	{
		double x = clamp((t - 150000.0) / (200000.0 - 150000.0), 0.0, 1.0);
		T=1.0-0.2*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (350000.0<=t and t<400000.0)
	{
		double x = clamp((t - 350000.0) / (400000.0 - 350000.0), 0.0, 1.0);
		T=0.8-0.2*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (550000.0<=t and t<600000.0)
	{
		double x = clamp((t - 550000.0) / (600000.0 - 550000.0), 0.0, 1.0);
		T=0.6-0.2*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (750000.0<=t and t<800000.0)
	{
		double x = clamp((t - 750000.0) / (800000.0 - 750000.0), 0.0, 1.0);
		T=0.4-0.2*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (950000.0<=t and t<1000000.0)
	{
		double x = clamp((t - 950000.0) / (1000000.0 - 950000.0), 0.0, 1.0);
		T=0.2-0.2*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	return T;
}

double Stepped(double t, double beta)
{
	double T=0.0;

	double Temp=1.0/beta;

	double interval=0.1-Temp;

	if (0.0<=t and t<=500000.0) {T=1.0;}
	else if (500000.0<t and t<=600000.0)
	{
		double x = clamp((t-500000.0) / (100000.0), 0.0, 1.0);
		T=1.0-0.25*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (600000.0<t and t<=1100000.0) {T=0.75;}
	else if (1100000.0<t and t<=1200000.0)
	{
		double x = clamp((t-1100000.0) / (100000.0), 0.0, 1.0);
		T=0.75-0.25*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (1200000.0<t and t<=1700000.0) {T=0.5;}
	else if (1700000.0<t and t<=1800000.0)
	{
		double x = clamp((t-1700000.0) / (100000.0), 0.0, 1.0);
		T=0.5-0.25*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (1800000.0<t and t<=2300000.0) {T=0.25;}
	else if (2300000.0<t and t<=2400000.0)
	{
		double x = clamp((t-2300000.0) / (100000.0), 0.0, 1.0);
		T=0.25-0.15*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (2400000.0<t and t<=3000000.0) {T=0.1;}
	else if (3000000.0<t and t<=4000000.0)
	{
		double x = clamp((t-3000000.0) / (1000000.0), 0.0, 1.0);
		T=0.1-interval*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (4000000.0<t and t<=5000000.0) {T=Temp;}
	
	return T;
}


double Critical_Temps(double t, double beta)
{
	double T=0.0,
	       Temp=1.0/beta;

	double interval=0.136-Temp;

	if (0.0<=t and t<=666667.0)
	{
		double x = clamp((t) / (666667.0), 0.0, 1.0);
		T=1.0-0.215*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (666667.0<t and t<=1666667.0) {T=0.785;}
	else if (1666667.0<=t and t<=2333334.0)
	{
		double x = clamp((t-1666667.0) / (2333334.0-1666667.0), 0.0, 1.0);
		T=0.785-0.649*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (2333334.0<t and t<=3333334.0) {T=0.136;}
	else if (3333334.0<=t and t<=4000000.0)
	{
		double x = clamp((t-3333334.0) / (4000000.0-3333334.0), 0.0, 1.0);
		T=0.136-interval*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (4000000.0<t and t<=5000000.0) {T=Temp;}
	
	return T;
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
		//crystal.rect_init(Num_Rep,Num_Rep,r);
		//crystal.ribbon_init(L,Num_Rep,r);
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

	ofstream MDVT;
	MDVT.open("mdvt.dat");

	ofstream AC;
	AC.open("AC.dat");

	lattice init;
	init = crystal;

	for (int t=restart_t; t<Time; t++)
	{
		// slope=10.0/(slp);
		// Temp=(1.0/cosh(w*slope*((double) t)))+Tf;
		//Temp=Step_Temp_Longer(t);
		Temp=Stepped(t,Tf);
		
		crystal.Metropolis(Temp,Edat,accepted, r);

		if (t==500000.0)
		{
			crystal.print_data("agg1.dat");

			stringstream sc_file;
			sc_file<<"sc_"<<t<<".dat";

			ofstream SC;
			SC.open(sc_file.str());
			crystal.spin_correlation(SC);
			SC.close();
		}

		if (t==1100000.0)
		{
			crystal.print_data("agg075.dat");
			stringstream sc_file;
			sc_file<<"sc_"<<t<<".dat";

			ofstream SC;
			SC.open(sc_file.str());
			crystal.spin_correlation(SC);
			SC.close();
		}
		if (t==1700000.0)
		{
			crystal.print_data("agg05.dat");
			stringstream sc_file;
			sc_file<<"sc_"<<t<<".dat";

			ofstream SC;
			SC.open(sc_file.str());
			crystal.spin_correlation(SC);
			SC.close();
		}
		if (t==2300000.0)
		{
			crystal.print_data("agg025.dat");
			stringstream sc_file;
			sc_file<<"sc_"<<t<<".dat";

			ofstream SC;
			SC.open(sc_file.str());
			crystal.spin_correlation(SC);
			SC.close();
		}
		if (t==3000000.0)
		{
			crystal.print_data("agg01.dat");
			stringstream sc_file;
			sc_file<<"sc_"<<t<<".dat";

			ofstream SC;
			SC.open(sc_file.str());
			crystal.spin_correlation(SC);
			SC.close();
		}

		// if (t%1000==0)
		// {
		// 	stringstream mov;
		// 	mov<<"mov_"<<t<<".dat";

		// 	crystal.print_data(mov.str());
		// 	rot_acc<<t<<" "<<accepted[1]/accepted[0]<<endl;
		// 	tran_acc<<t<<" "<<accepted[3]/accepted[2]<<endl;
		// 	glide_acc<<t<<" "<<accepted[5]/accepted[4]<<endl;
		// 	loc_acc<<t<<" "<<accepted[7]/accepted[6]<<endl;
		// }

		if (t%50000==0)
		{
			double autoC=0.0;
			for (int m=0; m<L*L; m++)
			{
				autoC+=((2*crystal.occupied(m)-1)*(2*init.occupied(m)-1));
			}
			AC<<t<<" "<<autoC<<endl;
		}

		//if (t%50000==0 and t>=1000000)
		if (t%50000==0)
		{
			stringstream inter;
			inter<<"Sys";
			crystal.print_data(inter.str());

			stringstream agg;
			agg<<"agg_"<<t;
			ofstream aggfile;
			aggfile.open(agg.str());
			HK clump(crystal);
			clump.Find_Cluster_periodic();

			double meanD=0.0;
			int lbl=clump.max_label();
			for (int i=1; i<=lbl; i++)
			{
				int ncl=clump.cluster_size(i);
				meanD += clump.distance_to_surface_periodic(i);
				if (ncl>=1)
				{
					aggfile<<ncl<<endl;
				}
			}
			aggfile.close();

			MDVT<<t<<" "<<meanD/N<<endl;
		}

		// if (t%50000==0 and t>=2000000)
		// {
		// 	stringstream sc_file;
		// 	sc_file<<"sc_"<<t<<".dat";

		// 	ofstream SC;
		// 	SC.open(sc_file.str());
		// 	crystal.spin_correlation(SC);
		// 	SC.close();
		// }
	}

	MDVT.close();
	AC.close();

	ofstream SC;
	SC.open("sc_final.dat");

	crystal.spin_correlation(SC);
	
	SC.close();

	// for (int t=0; t<=end_time; t++)
	// {
	// 	crystal.Metropolis(0.0,Edat,accepted, r);
	// }

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

	HK Aggregate(crystal);
	Aggregate.Find_Cluster_periodic();
	Aggregate.print_cluster();
	int clusters=Aggregate.cluster_count();

	for (int n=1; n<=Aggregate.max_label();n++)
	{
		int size = Aggregate.cluster_size(n);
		if (size == 0) {continue;}
		vector<double> md = Aggregate.mean_distance_to_surface_periodic(n);
		inf<<"	Mean Distance to Surface: "<<md[0]<<" STD: "<<md[1]<<endl;
	}

	inf<<"Clusters: "<<clusters<<endl;

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
	inf<<"J="<<J<<" Sigma="<<Sigma<<" K="<<K<<" f="<<f<<endl;
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