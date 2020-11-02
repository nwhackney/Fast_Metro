#include "../include/lattice.hpp"
#include <stdlib.h>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <limits>

double Box_Muller(double mu, double sigma, gsl_rng * rng)
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
	   //u1 = rand() * (1.0 / RAND_MAX);
	   //u2 = rand() * (1.0 / RAND_MAX);
	 	u1 = gsl_rng_uniform(rng);
	 	u2 = gsl_rng_uniform(rng);
	 }
	while ( u1 <= epsilon );

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

void lattice::Metropolis(double T, std::ofstream &Efile, std::vector<double> &accepted, gsl_rng * rng) // Accepted: (move 0 tried, move 0 accepted, move 1 tried, move 1 accepted ...)
{
	double Total_Energy=H();
	Efile<<Total_Energy<<std::endl;
	for (int i=0; i<N; i++)
	{
		int n=occ[i];
		double E=H_local(n);

		std::vector<site> saved = lattice;

		//int flag=rand()%4;
		int flag=gsl_rng_get(rng)%4;
		//int flag=0;

		if (flag==0) // Rotation
		{
			accepted[0]+=1.0;
			
			//double width = 3.1415962;
			//double delta=2.0*(gsl_rng_uniform(rng)-0.5);
			//double theta=lattice[n].angle+delta*width;

			double width = 0.1*exp(1.5*T);
			double theta = Box_Muller(lattice[n].angle,width,rng); // Maybe should switch to a gsl gaussian random distribution, instead of "hacking" the Box_Muller Function (twould be more elegante)
			rotate(n,theta);
			double Trial_E=H_local(n);
			double delE = Trial_E - E;
			
			double alpha = gsl_rng_uniform(rng);
			double U= exp(-1*delE/T);
			if (alpha > fmin(1.0,U))
			{
				lattice=saved;
			}
			else {accepted[1]+=1.0;}
		}
		else if (flag==1) // Translation
		{
			accepted[2]+=1.0;

			int unocc= gsl_rng_get(rng)%(V-N);
			int m=vac[unocc];

			double phi=lattice[n].angle;

			site Null; Null.occ=0; Null.angle=0.0;
			site Spin; Spin.occ=1; Spin.angle=phi;

			lattice.at(n)=Null;
			lattice.at(m)=Spin;

			double Trial_E=H_local(m);

			double delE = Trial_E - E;
			//double alpha = ((double) rand()/(double) RAND_MAX);
			double alpha = gsl_rng_uniform(rng);

			double U= exp(-1.0*delE/T);
			if (alpha > fmin(1.0,U))
			{
				lattice.clear();
				lattice=saved;
			}
			else
			{
				occ[i]=m;
				vac[unocc]=n;
				accepted[3]+=1.0;
			}

		}
		else if (flag==2) // Translation + Rotation
		{
			accepted[4]+=1.0;
			
			int unocc= gsl_rng_get(rng)%(V-N);
			int m=vac[unocc];
			
			double theta = gsl_rng_uniform(rng)*6.283185307179586;
			
			site Null; Null.occ=0; Null.angle=0.0;
			site Spin; Spin.occ=1; Spin.angle=theta;

			lattice.at(n)=Null;
			lattice.at(m)=Spin;
			
			double Trial_E=H_local(m);
			
			double delE = Trial_E - E;
			double alpha = gsl_rng_uniform(rng);
			
			double U= exp(-1*delE/T);
			if (alpha > fmin(1.0,U))
			{
				lattice=saved;
			}
			else
			{
				occ[i]=m;
				vac[unocc]=n;
				accepted[5]+=1.0;
			}

		}
		else // Local Translation
		{
			int slide=gsl_rng_get(rng)%4;
			int slot;
			
			if (slide==0)                     // Up
			{
				slot=(n-L)%V;
				if (slot<0){slot=slot+V;}
			}
			else if (slide==1){slot=(n+L)%V;} // Down
			else if (slide==2)                // Left
			{
				slot=(n-1)%V;
				if (slot<0){slot=slot+V;}
			}
			else              {slot=(n+1)%V;} // Right
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
			
			accepted[6]+=1.0;
			
			std::vector<int> occ_saved = occ;
			std::vector<int> vac_saved = vac;
			
			double theta = lattice[n].angle;
			
			site Null; Null.occ=0; Null.angle=0.0;
			site Spin; Spin.occ=1; Spin.angle=theta;
			
			lattice.at(n)=Null;
			lattice.at(slot)=Spin;
			
			occ[i]=slot; vac[m]=n;
			
			double Trial_E=H_local(slot);
			
			double delE = Trial_E - E;
			//double alpha = ((double) rand()/(double)RAND_MAX);
			double alpha=gsl_rng_uniform(rng);

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


void lattice::Spin_Metropolis(double T, std::ofstream &Efile, std::vector<double> &accepted, gsl_rng * rng) // Accepted: (move 0 tried, move 0 accepted, move 1 tried, move 1 accepted ...)
{
	double Total_Energy=H();
	Efile<<Total_Energy<<std::endl;
	for (int i=0; i<N; i++)
	{
		int n=occ[i];
		double E=H_local(n);

		std::vector<site> saved = lattice;

		//int flag=rand()%4;
		
		accepted[0]+=1.0;
		
		//double width = 3.1415962;
		//double delta=2.0*(gsl_rng_uniform(rng)-0.5);
		//double theta=lattice[n].angle+delta*width;

		double width = 0.1*exp(1.5*T);
		double theta = Box_Muller(lattice[n].angle,width,rng); // Maybe should switch to a gsl gaussian random distribution, instead of "hacking" the Box_Muller Function (twould be more elegante)
		rotate(n,theta);
		double Trial_E=H_local(n);
		double delE = Trial_E - E;
		
		double alpha = gsl_rng_uniform(rng);
		double U= exp(-1*delE/T);
		if (alpha > fmin(1.0,U))
		{
			lattice=saved;
		}
		else {accepted[1]+=1.0;}
		
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


void lattice::init(int lattice_size, int Number, gsl_rng * rng)
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

	// for (int i=0; i<V; i++)
	// {
	// 	double theta = gsl_rng_uniform(rng)*6.283185307179586;

	// 	lattice[n].occ=1;
	// 	lattice[n].angle=theta;
	// }

	int count=0;
	while(count!=N)
	{
		int n=gsl_rng_get(rng)%V;
		if (n<0){n=n+V;}
		if (lattice[n].occ==1) {continue;}

		double theta = gsl_rng_uniform(rng)*6.283185307179586;

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

void lattice::rect_init(int Lx, int Ly, gsl_rng * rng)
{
	
	if (Lx>=Ly)
	{
		V=(Lx+2)*(Lx+2);
		L=Lx+2;
	}
	else if (Ly>Lx)
	{
		V=(Ly+2)*(Ly*2);
		L=Ly+2;
	}

	N=Lx*Ly;

	site Null;
	Null.occ=0;
	Null.angle=0.0;

	lattice.resize(V,Null);
	occ.resize(N);
	vac.resize(V-N);

	int count=0;
	for (int i=0; i<Lx; i++)
	{
		for (int j=0; j<Ly; j++)
		{
			int n = (i+1)*L+(j+1);
			double theta = gsl_rng_uniform(rng)*6.283185307179586;

			lattice[n].occ=1;
			lattice[n].angle=theta;

			occ[count]=n;
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
			//out << x << " " << y << " " << lattice[i].angle << std::endl;
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
	out<<"set xrange [0:75]"<<std::endl;
	out<<"set yrange [0:75]"<<std::endl;
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
	//f=0.5*frustration; // Factor of 1/2 to match continuum calculation, where I dropped the 2pi factor
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
	int up=(i-L)%V,
	    down=(i+L)%V,
	    left=(i-1)%V,
	    right=(i+1)%V;

    if (up<0){up=up+V;}
    if (left<0){left=left+V;}

    double x = (double) (i%L);
    double y = (double) (i/L);

    double weight_up=((double) lattice[i].occ)*((double) lattice[up].occ);
    double H1 = (-1.0*J*cos(lattice[i].angle-lattice[up].angle+f*x)-K)*weight_up;                    			    // Up Bond

    double weight_down=((double) lattice[i].occ)*((double) lattice[down].occ);
    double H2 = (-1.0*J*cos(lattice[i].angle-lattice[down].angle-f*x)-K)*weight_down; 		                        // Down Bond

    double weight_left=((double) lattice[i].occ)*((double) lattice[left].occ);
    double H3 =(-1.0*J*cos(lattice[i].angle-lattice[left].angle-f*y)-K)*weight_left;	                                  // Left Bond

    double weight_right=((double) lattice[i].occ)*((double) lattice[right].occ);
    double H4 =(-1.0*J*cos(lattice[i].angle-lattice[right].angle+f*y)-K)*weight_right;                    	         // Right Bond

	double H = H1 + H2 + H3 + H4;

	return H;
}

double lattice::H()
{
	double H=0.0;

	for (int i=0; i<V; i++)
	{
		int right=(i+1)%V,
		    down=(i+L)%V;

		    H+=(-1.0*J*cos(lattice[i].angle-lattice[right].angle+f*(i/L))-K)*lattice[i].occ*lattice[right].occ;
		    H+=(-1.0*J*cos(lattice[i].angle-lattice[down].angle-f*(i%L))-K)*lattice[i].occ*lattice[down].occ;
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

    if (up<0){up=up+V;}
    if (left<0){left=left+V;}

    H+=1.0*J*cos(lattice[i].angle-lattice[up].angle+f*(i%L))*lattice[i].occ*lattice[up].occ; 			    // Up Bond
    H+=-1.0*J*cos(lattice[i].angle-lattice[down].angle-f*(i%L))*lattice[i].occ*lattice[down].occ; 		    // Down Bond
    H+=-1.0*J*cos(lattice[i].angle-lattice[left].angle-f*(i/L))*lattice[i].occ*lattice[left].occ;		    // Left Bond
    H+=-1.0*J*cos(lattice[i].angle-lattice[right].angle+f*(i/L))*lattice[i].occ*lattice[right].occ; 	    // Right Bond

    int neighbor=lattice[up].occ+lattice[down].occ+lattice[left].occ+lattice[right].occ;
    double num= (double) neighbor;


	return H/num;
}
