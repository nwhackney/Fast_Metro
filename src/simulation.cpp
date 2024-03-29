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
	const toml::Value* ACF = v.find("AutoCorrelation");
	
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
	acf=ACF->as<string>();
}

//Temperature Schedules Below

float clamp(float x, float lowerlimit, float upperlimit)
{
  if (x < lowerlimit)
    x = lowerlimit;
  if (x > upperlimit)
    x = upperlimit;
  return x;
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
	else if (4000000.0<t and t<=6000000.0) {T=Temp;}
	
	return T;
}

double Stepped_Long(double t, double beta)
{
	double T=0.0;

	double Temp=1.0/beta;

	double interval=0.1-Temp;

	if (0.0<=t and t<=1000000.0) {T=1.0;}
	else if (1000000.0<t and t<=1200000.0)
	{
		double x = clamp((t-1000000.0) / (200000.0), 0.0, 1.0);
		T=1.0-0.25*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (1200000.0<t and t<=2200000.0) {T=0.75;}
	else if (2200000.0<t and t<=2400000.0)
	{
		double x = clamp((t-2200000.0) / (200000.0), 0.0, 1.0);
		T=0.75-0.25*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (2400000.0<t and t<=3400000.0) {T=0.5;}
	else if (3400000.0<t and t<=3600000.0)
	{
		double x = clamp((t-3400000.0) / (200000.0), 0.0, 1.0);
		T=0.5-0.25*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (3600000.0<t and t<=4600000.0) {T=0.25;}
	else if (4600000.0<t and t<=4800000.0)
	{
		double x = clamp((t-4600000.0) / (200000.0), 0.0, 1.0);
		T=0.25-0.15*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (4800000.0<t and t<=5800000.0) {T=0.1;}
	else if (5800000.0<t and t<=8000000.0)
	{
		double x = clamp((t-5800000.0) / (2200000.0), 0.0, 1.0);
		T=0.1-interval*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (8000000.0<t and t<=12000000.0) {T=Temp;}
	
	return T;
}

double Elongated(double t, double beta)
{
	double T=0.0;

		double Temp=1.0/beta;

		double interval=(0.1-Temp)/3.0;

		double TS1=0.1-interval,
			  TS2=TS1-interval,
			  TS3=TS2-interval;

		if (0.0<=t and t<=250000.0) {T=1.0;}
		else if (250000.0<t and t<=260000.0)
		{
			double x = clamp((t-250000.0) / (100000.0), 0.0, 1.0);
			T=1.0-0.25*x*x*x*(x*(x*6.0-15.0)+10.0);
		}
		else if (260000.0<t and t<=510000.0) {T=0.75;}
		else if (510000.0<t and t<=520000.0)
		{
			double x = clamp((t-510000.0) / (100000.0), 0.0, 1.0);
			T=0.75-0.25*x*x*x*(x*(x*6.0-15.0)+10.0);
		}
		else if (520000.0<t and t<=770000.0) {T=0.5;}
		else if (770000.0<t and t<=780000.0)
		{
			double x = clamp((t-770000.0) / (100000.0), 0.0, 1.0);
			T=0.5-0.25*x*x*x*(x*(x*6.0-15.0)+10.0);
		}
		else if (780000.0<t and t<=1030000.0) {T=0.25;}
		else if (1030000.0<t and t<=1040000.0)
		{
			double x = clamp((t-1030000.0) / (100000.0), 0.0, 1.0);
			T=0.25-0.15*x*x*x*(x*(x*6.0-15.0)+10.0);
		}
		else if (1040000.0<t and t<=1500000.0) {T=0.1;}
		else if (1500000.0<t and t<=2000000.0)
		{
			double x = clamp((t-1500000.0) / (500000.0), 0.0, 1.0);
			T=0.1-interval*x*x*x*(x*(x*6.0-15.0)+10.0);
		}
		else if (2000000.0<t and t<=3000000.0) {T=TS1;}
		else if (3000000.0<t and t<=3500000.0)
		{
			double x = clamp((t-3000000.0) / (500000.0), 0.0, 1.0);
			T=TS1-interval*x*x*x*(x*(x*6.0-15.0)+10.0);
		}
		else if (3500000.0<t and t<=4500000.0) {T=TS2;}
		else if (4500000.0<t and t<=5000000.0)
		{
			double x = clamp((t-4500000.0) / (500000.0), 0.0, 1.0);
			T=TS2-interval*x*x*x*(x*(x*6.0-15.0)+10.0);
		}
		else if (5000000.0<t and t<=8000000.0) {T=TS3;}
	
	return T;
}

double No_Step(double t, double beta)
{
	double T=0.0,
	       Temp=1.0/beta;

	double interval=1.0;

	if (0.0<=t and t<=2000000.0)
	{
		double x = clamp((t) / (2000000.0), 0.0, 1.0);
		T=1.0-interval*x*x*x*(x*(x*6.0-15.0)+10.0);
	}
	else if (2000000.0<t and t<=3000000) {T=0.0;}
	
	return T;
}

// Temperature Schedules Above

// Auto-Correlation Functions ~Attempts~ Below

void Cluster_AC(lattice &tc, int t, vector<double> &time, vector<double> &avg)
{
	HK current(tc);
	current.Find_Cluster_periodic();

	double current_avg=0.0;
	int current_label=current.max_label();
	double current_num=(double) current.cluster_count();

	for (int i=1; i<=current_label;i++)
	{
		int cnc=current.cluster_size(i);
		current_avg+=(double) cnc;
	}
	current_avg=current_avg/current_num;

	time.push_back(t); avg.push_back(current_avg);	
}

void Uni_AC(lattice &tc, vector< vector<double> > &avg)
{
	HK current(tc);
	current.Find_Cluster_periodic();

	for(int n=0; n<tc.how_many(); n++)
	{
		int index=tc.occ_to_index(n);
		int label=current.index_to_label(index);

		int size=current.cluster_size(label);
		//cout<<index<<" "<<tc.occupied(index)<<" "<<label<<" "<<size<<endl;

		avg[n].push_back(size);	
	}
}

// Auto-Correlation Functions ~Attempts~ Above

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
		//crystal.init(L,N,r);
		crystal.square_init(L,N,Num_Rep,r);
		//crystal.ribbon_init(L,Num_Rep,r);
		//crystal.hex_init(L,Num_Rep,r);
		crystal.print_data("init");	
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

	ofstream Navg;
	Navg.open("Navg.dat");

	lattice init;
	init = crystal;

	// Temp Looking at Agg of initialization
	ofstream Taggfile;
	Taggfile.open("agg_init.dat");
	HK clump(crystal);
	clump.Find_Cluster_periodic();

	double MDD=0.0;
	int lab=clump.max_label();
	for (int i=1; i<=lab; i++)
	{
		int numc=clump.cluster_size(i);
		MDD += clump.distance_to_surface_periodic(i);
		if (numc>=1)
		{
			Taggfile<<numc<<endl;
		}
	}
	Taggfile.close();
	// End Temp

	vector<double> slice;
	// vector< vector<double> > agav;
	// agav.resize(N);
	
	//***************** Below commented because I don't want to deal
	// FILE *h_fp = fopen("cluster_hist.out", "w");
	// FILE *d_fp = fopen("cluster_dat.out", "w");
	// int global_hist[N] = {0};
	//****************
	
	for (int t=restart_t; t<Time; t++)
	{
		// Selects annealing schedule and calculates temperature
		// slope=10.0/(slp);
		// Temp=(1.0/cosh(w*slope*((double) t)))+Tf;
		//Temp=Step_Temp_Longer(t);
		//Temp=QUENCH(t,Tf);
		Temp=No_Step(t,Tf);
		//Temp=1.0/Tf;
		
		crystal.Spin_Metropolis(Temp,Edat,accepted, r); // runs monte carlo step
		//Uni_AC(crystal, agav);

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

		if (t > 100000 && t%50000 == 0)
		{

			//************************** Below Commented because I don't want to deal with it

			// stringstream inter;
			// inter<<"Sys";
			// crystal.print_data(inter.str());

			// stringstream agg;
			// agg<<"agg_"<<t;
			// crystal.print_data(agg.str());


			// // Commented For Movie

			// ofstream aggfile;
			// aggfile.open(agg.str());
			// HK clump(crystal);				// Hoshen-Kopelman (HK) object for finding/counting clusters
			// clump.Find_Cluster_periodic();	// Identifies clusters

			// double avg_size = 0.0, avg_size_2 = 0.0;
			// // max_label = num_clusters?
			// int max_label = clump.max_label();
			// int num_clusters = clump.cluster_count();
			// for (int j = 1; j <= max_label; j++)		// Finds average cluster size
			// {
			// 	int size = clump.cluster_size(j);
			// 	avg_size += size;
			// 	avg_size_2 += size *size;
			// }
			// avg_size /= num_clusters;
			// avg_size_2 /= num_clusters;
			// Navg << t << " " << avg_size << endl;
			// fprintf(d_fp, "%15.15e %15.15e %d\n", avg_size, avg_size_2, num_clusters);
			// fflush(d_fp);

			// double meanD=0.0;
			// int lbl=clump.max_label();
			// int cluster_size_hist[N] = {0};
			// for (int i=1; i<=lbl; i++)		// Finds Mean distance to surface
			// {
			// 	int ncl=clump.cluster_size(i);
			// 	double dts = clump.distance_to_surface_periodic(i);
			// 	meanD += dts;
			// 	if (ncl>=1)
			// 	{
			// 		aggfile<<ncl<<" "<<dts/((double) ncl)<<endl;
			// 		cluster_size_hist[ncl]++;
			// 		global_hist[ncl]++;
			// 	}
			// }

			// for (int i = 1; i <= N; i++) {
			// 	if (cluster_size_hist[i] > 0)
			// 		fprintf(h_fp, "%d ", cluster_size_hist[i]);
			// }
			// fprintf(h_fp, "\n");

			// for (int i = 1; i <= N; i++) {
			// 	if (cluster_size_hist[i] > 0)
			// 		fprintf(h_fp, "%d ", i);
			// }
			// fprintf(h_fp, "\n");
			// fflush(h_fp);
			// aggfile.close();

			// MDVT<<t<<" "<<meanD/N<<endl;

			// Above Commented for Movie
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

		//*********************
	}

     //*********************** Below Commented because I don't want to deal with is
	
	// for (int i = 1; i <= N; i++) {
	// 	if (global_hist[i] > 0)
	// 		fprintf(h_fp, "%d ", global_hist[i]);
	// }
	// fprintf(h_fp, "\n");

	// for (int i = 1; i <= N; i++) {
	// 	if (global_hist[i] > 0)
	// 		fprintf(h_fp, "%d ", i);
	// }
	// fclose(d_fp);
	// fclose(h_fp);

	//*****************************
	
	// ofstream AutoC;
	// AutoC.open("Auto.dat");
	// double chi_zero=0.0;
	// for (int i=0; i<=Time; i++)
	// {
	// 	int diff=Time-i;
		
	// 	double coeff=1.0/((double) diff);
	// 	double chi=0.0,
	//             nChi=0.0;

	//      //int l=end_time;
	// 	for (int l=0; l<N; l++)
	// 	{
	// 		for (int j=0; j<=diff; j++)
	// 		{
	// 			chi+=coeff*agav[l][j]*agav[l][j+i];
	// 			for (int k=0; k<=diff; k++)
	// 			{
	// 				chi-=coeff*coeff*agav[l][j]*agav[l][k+i];
	// 			}
	// 		}
	// 		if(i==0){chi_zero=chi;}
	// 		nChi+=chi/chi_zero;
	// 	}
	// 	AutoC<<i<<" "<<nChi/N<<endl;
	// }
	// AutoC.close();

	//MDVT.close();
	Navg.close();

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

	ofstream inf;						// Prints Info File
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
	Aggregate.clusters_labelled();
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
