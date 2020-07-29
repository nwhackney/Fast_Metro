#include "../include/lattice.hpp"
#include <stdlib.h>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <limits>

double Box_Muller(double mu, double sigma)
{
	static const double epsilon = std::numeric_limits<double>::min();
	static const double two_pi = 2.0*3.1416;

	thread_local double z1;
	thread_local bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	double u1, u2;
	do
	 {
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	 }
	while ( u1 <= epsilon );

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

void lattice::Metropolis(double T, std::ofstream &Efile, std::vector<double> &accepted) // Accepted: (move 0 tried, move 0 accepted, move 1 tried, move 1 accepted ...)
{
	double Total_Energy=H();
	Efile<<Total_Energy<<std::endl;

	for (int i=0; i<N; i++)
	{
		int n=occ[i];
		double E=H_local(n),
		       Trial_E;

		std::vector<site> saved = lattice;

		int flag=rand()%1;

		if (flag==0) // Rotation
		{
			accepted[0]+=1.0;

			double width = 0.4*exp(T);
			double theta = Box_Muller(lattice[n].angle,width);

			rotate(n,theta);
			Trial_E=H_local(n);

			double delE = Trial_E - E;
			double alpha = ((double) rand()/(double)RAND_MAX);

			double U= exp(-1*delE/T);
			if (alpha > fmin(1.0,U))
			{
				lattice=saved;
			}
			else {accepted[1]+=1.0;}
		}
		// else if (flag==1) // Translation
		// {
		// 	accepted[2]+=1.0;

		// 	int unocc= rand()%(V-N);
		// 	int m=vac[unocc];

		// 	double phi=lattice[n].angle;
		// 	flip(n);
		// 	rotate(n,0.0);
		// 	flip(m);
		// 	rotate(m,phi);

		// 	Trial_E=H_local(m);

		// 	double delE = Trial_E - E;
		// 	double alpha = ((double) rand()/(double)RAND_MAX);

		// 	double U= exp(-1*delE/T);
		// 	if (alpha > fmin(1.0,U))
		// 	{
		// 		lattice=saved;
		// 	}
		// 	else
		// 	{
		// 		occ[i]=m;
		// 		vac[unocc]=n;
		// 		accepted[3]+=1.0;
		// 	}

		// }
		// else if (flag==2) // Translation + Rotation
		// {
		// 	accepted[4]+=1.0;

		// 	int unocc= rand()%(V-N);
		// 	int m=vac[unocc];

		// 	flip(n);
		// 	rotate(n,0.0);
		// 	flip(m);

		// 	double theta = ((double) rand()*(6.28)/(double)RAND_MAX);
		// 	rotate(m,theta);

		// 	Trial_E=H_local(m);

		// 	double delE = Trial_E - E;
		// 	double alpha = ((double) rand()/(double)RAND_MAX);

		// 	double U= exp(-1*delE/T);
		// 	if (alpha > fmin(1.0,U))
		// 	{
		// 		lattice=saved;
		// 	}
		// 	else
		// 	{
		// 		occ[i]=m;
		// 		vac[unocc]=n;
		// 		accepted[5]+=1.0;
		// 	}

		//}
		else if (flag==1) // Local Translation
		{
			accepted[6]+=1.0;

			int slide=rand()%4;
			int slot;

			if (slide==0){slot=(n-L)%V;}      // Up
			else if (slide==1){slot=(n+L)%V;} // Down
			else if (slide==2){slot=(n-1)%V;} // Left
			else if (slide==3){slot=(n+1)%V;} // Right

			int VS=V-N;
			int m=-1;
			for (int j=0; j<VS; j++)
			{
				if (vac[j]==slot)
				{
					m=j;
				}
			}

			if (m==-1) {continue;}

			std::vector<int> occ_saved = occ;
			std::vector<int> vac_saved = vac;

			double theta = lattice[n].angle;
			lattice[n].angle=0.0; lattice[n].occ=0;
			lattice[slot].angle=theta; lattice[slot].occ=1;

			occ[i]=slot; vac[m]=n;

			Trial_E=H_local(slot);

			double delE = Trial_E - E;
			double alpha = ((double) rand()/(double)RAND_MAX);

			double U= exp(-1*delE/T);
			if (alpha > fmin(1.0,U))
			{
				lattice=saved;
				occ=occ_saved;
				vac=vac_saved;
			}
			else {accepted[7]+=1.0;}
		}

	}
}

void lattice::check()
{
	//std::cout<<"Volume: "<<V<<" Array Size: "<<lattice.size()<<std::endl;
	int c=0;
	for (int i=0; i<N; i++)
	{
		int test=occ[i];
		if (lattice[test].occ==1)
		{
			c++;
		}
	}
	std::cout<<c<<std::endl;
}


void lattice::init(int lattice_size, int Number)
{

	L=lattice_size;
	N=Number;

	V=L*L;

	site Null;
	Null.occ=0;
	Null.angle=0.0;

	lattice.resize(V,Null);
	occ.resize(N);
	vac.resize(V-N);

	int count=0;
	while(count!=N)
	{
		int n=(rand()+getpid())%V;
		if (lattice[n].occ==1) {continue;}

		double theta = ((double) rand()*(6.28)/(double)RAND_MAX);

		lattice[n].occ=1;
		lattice[n].angle=theta;

		occ[count]=n;
		count++;
	}

	int m=0;
	for (int i=0; i<V; i++)
	{
		if (lattice[i].occ==0)
		{
			vac[m]=i;
			m++;
		}
	}
	//std::cout<<occ.size()<<" "<<vac.size()<<" "<<occ.size()+vac.size()<<std::endl;
}

void lattice::restart(int Length, int Number, std::string infile)
{

	L=Length;
	N=Number;
	V=L*L;
	
	std::ifstream in;
	in.open(infile);

	if (!in)
	{
     	std::cout << "Unable to open file";
     	exit(1); // terminate with error
     }

     std::vector<double> raw;
     std::string x;

     while (in>>x)
     {
     	std::stringstream item;
     	item<<x;
     	double i=0;
     	item >> i;
     	raw.push_back(i);
     }

     in.close();

     site Null;
	Null.occ=0;
	Null.angle=0.0;

	lattice.resize(V,Null);
	occ.resize(N);
	vac.resize(V-N);

	int raw_size=raw.size();

	int count=0;

     for (int i=0; i<=raw_size-4; i+=4)
     {
          if (raw[i+2]!=-1.0)
          {
          	int n=(int) raw[i];
          	int m=(int) raw[i+1];
          	int pos=n*L+m;
          	double theta=(double) raw[i+2];
               lattice[pos].occ=1;
               lattice[pos].angle=theta;

               occ[count]=pos;
               count++;
          }
     }

     int m=0;
	for (int i=0; i<V; i++)
	{
		if (lattice[i].occ==0)
		{
			vac[m]=i;
			m++;
		}
	}
}

void lattice::print_data(std::string file_name)
{
	std::stringstream file;
	file<<file_name<<"_data.dat";

	std::ofstream out;
	out.open(file.str());

	for (int i=0; i<V; i++)
	{
		int x=i%L;
		int y=i/L;

		if (lattice[i].occ==0)
		{
			out<<x<<" "<<y<<" "<<" -1 "<<" 0"<<std::endl;
		}
		else
		{
			out << x << " " << y << " " << lattice[i].angle << " " << strain(i) << std::endl;
		}
	}
}

void lattice::print_gnu(std::string file_name)
{

	std::stringstream file;
	file<<file_name<<".p";

	std::ofstream out;
	out.open(file.str());

	std::stringstream png;
	png<<file_name<<".png";

	out<<"set terminal png"<<std::endl;
	out<<"set output '"<<png.str()<<"'"<<std::endl;
	out<<"set key off"<<std::endl;
	out<<"set xrange [0:27]"<<std::endl;
	out<<"set yrange [0:27]"<<std::endl;
	out<<"set style arrow 1 head filled size screen 0.03,15 ls 2"<<std::endl;

	double d=2.5;
	double theta,
	       x,
	       y,
	       dx,
	       dy;

	for (int i=0; i<V; i++)
	{
			if (lattice[i].occ==1)
			{
				int n=i%L,
				    m=i/L;

				x=(n+1)*d;
				y=(m+1)*d;

				theta=lattice[i].angle;
				dx=cos(theta);
				dy=sin(theta);

				out<<"set arrow from "<<x<<","<<y<<" to "<<x+dx<<","<<y+dy<<" as 1 lc 'black'"<<std::endl;
			}
	}

	out<<"plot NaN"<<std::endl;
}

void lattice::set_const(double j, double sig, double frustration)
{
	J=j;
	K=sig;
	f=frustration*3.1415962; //Why is this multiplied by Pi? Need to make sure that is correct...
}

void lattice::rotate(int i, double theta)
{
	lattice[i].angle=theta;
}

void lattice::flip(int i)
{
	if (lattice[i].occ==1){lattice[i].occ=0;}
	else {lattice[i].occ=1;}
}

int lattice::occupied(int i)
{
	return lattice[i].occ;
}

int lattice::how_many()
{
	return N;
}

double lattice::angle(int i)
{
	return lattice[i].angle;
}

double lattice::H_local(int i)
{
	double H=0.0;

	int up=(i-L)%V,
	    down=(i+L)%V,
	    left=(i-1)%V,
	    right=(i+1)%V;

	    H+=(J*cos(lattice[i].angle-lattice[up].angle+f*(i%L))+K)*lattice[i].occ*lattice[up].occ; 			    // Up Bond
	    H+=(J*cos(lattice[i].angle-lattice[down].angle-f*(i%L))+K)*lattice[i].occ*lattice[down].occ; 		    // Down Bond
	    H+=(J*cos(lattice[i].angle-lattice[left].angle-f*(i/L))+K)*lattice[i].occ*lattice[left].occ;		    // Left Bond
	    H+=(J*cos(lattice[i].angle-lattice[right].angle+f*(i/L))+K)*lattice[i].occ*lattice[right].occ; 	         // Right Bond


	return H;
}

double lattice::H()
{
	double H=0.0;

	for (int i=0; i<V; i++)
	{
		int right=(i+1)%V,
		    down=(i+L)%V;

		    H+=(J*cos(lattice[i].angle-lattice[right].angle+f*(i/L))+K)*lattice[i].occ*lattice[right].occ;
		    H+=(J*cos(lattice[i].angle-lattice[down].angle-f*(i%L))+K)*lattice[i].occ*lattice[down].occ;
	}

	return H;
}

double lattice::strain(int i)
{
	double H=0.0;

	int up=(i-L)%V,
	    down=(i+L)%V,
	    left=(i-1)%V,
	    right=(i+1)%V;

	    H+=J*cos(lattice[i].angle-lattice[up].angle+f*(i%L))*lattice[i].occ*lattice[up].occ; 			    // Up Bond
	    H+=J*cos(lattice[i].angle-lattice[down].angle-f*(i%L))*lattice[i].occ*lattice[down].occ; 		    // Down Bond
	    H+=J*cos(lattice[i].angle-lattice[left].angle-f*(i/L))*lattice[i].occ*lattice[left].occ;		    // Left Bond
	    H+=J*cos(lattice[i].angle-lattice[right].angle+f*(i/L))*lattice[i].occ*lattice[right].occ; 	         // Right Bond

	    int neighbor=lattice[up].occ+lattice[down].occ+lattice[left].occ+lattice[right].occ;
	    double num= (double) neighbor;


	return H/num;
}
