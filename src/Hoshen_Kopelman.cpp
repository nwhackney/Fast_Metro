#include "../include/lattice.hpp"
#include "../include/Hoshen_Kopelman.hpp"
//#include "Skeleton.cpp"
#include <fstream>
#include <vector>

using namespace std;

int find(int x, std::vector<int> &labels)
{
	int y = x;
	while (labels[y] != y)
		y = labels[y];

	while (labels[x] != x)
	{
		int z = labels[x];
		labels[x] = y;
		x = z;
	}
	return y;
}

int unionize(int x, int y, std::vector<int> &labels)
{
	return labels[find(x, labels)] = find(y, labels);
}

HK::HK(lattice init)
{
	system=init;
	L=system.size();
	N_Spins=system.how_many();	
}

void HK::Find_Cluster()
{
	for (int i=0; i<N_Spins; i++)
	{
		labels.push_back(i);
	}

	matrix.resize(L);
	for (int d=0; d<L; d++)
	{
		matrix[d].resize(L,0);
	}

	for (int n=0; n<L; n++)
	{
		for (int m=0; m<L; m++)
		{
			if (system.occupied(n*L+m) == 0) {continue;}
			else {matrix[n][m] = 1;}
		}
	}

	for (int i=0; i<L; i++)
	{
		for (int j=0; j<L; j++)
		{
			if (matrix[i][j])
			{                        // if occupied ...
				int up = (i==0 ? 0 : matrix[i-1][j]);    //  look up  
				int left = (j==0 ? 0 : matrix[i][j-1]);  //  look left
				
				switch (!!up + !!left)
				{
		  			case 0:
		  				labels[0] ++;
						assert(labels[0] < N_Spins);
						labels[labels[0]] = labels[0];
		  				matrix[i][j] = labels[0];      // a new cluster
		  				break;
		  
					case 1:                              // part of an existing cluster
						matrix[i][j] = max(up,left);    // whichever is nonzero is labelled
						break;
		  
					case 2:                              // this site binds two clusters
						matrix[i][j] = unionize(up, left, labels);
						break;
				}
			}
		}
	}

	for (int i=0; i<L; i++)
	{
		for (int j=0; j<L; j++)
		{
			if (matrix[i][j]==0) {continue;}
			matrix[i][j]=find(matrix[i][j],labels);
		}
	}
}

void HK::Find_Cluster_periodic()
{
	for (int i=0; i<N_Spins; i++)
	{
		labels.push_back(i);
	}

	matrix.resize(L);
	for (int d=0; d<L; d++)
	{
		matrix[d].resize(L,0);
	}

	for (int n=0; n<L; n++)
	{
		for (int m=0; m<L; m++)
		{
			if (system.occupied(n*L+m) == 0) {continue;}
			else {matrix[n][m] = 1;}
		}
	}

	for (int i=0; i<L; i++)
	{
		for (int j=0; j<L; j++)
		{
			if (matrix[i][j])
			{                        // if occupied ...
				int up = (i==0 ? matrix[L-1][j] : matrix[i-1][j]);    //  look up  
				int left = (j==0 ? matrix[i][L-1] : matrix[i][j-1]);  //  look left
				
				switch (!!up + !!left)
				{
		  			case 0:
		  				labels[0]++;
						assert(labels[0] < N_Spins);
						labels[labels[0]] = labels[0];
		  				matrix[i][j] = labels[0];      // a new cluster
		  				break;
		  
					case 1:                              // part of an existing cluster
						matrix[i][j] = max(up,left);    // whichever is nonzero is labelled
						break;
		  
					case 2:                              // this site binds two clusters
						matrix[i][j] = unionize(up, left, labels);
						break;
				}

				for (int i=0; i<L; i++)
				{
					if (matrix[i][0]!=0 and matrix[i][L-1]!=0) {matrix[i][L-1]=unionize(matrix[i][0],matrix[i][L-1],labels);}
				}

				for (int j=0; j<L; j++)
				{
					if (matrix[0][j]!=0 and matrix[L-1][j]!=0) {matrix[L-1][j]=unionize(matrix[0][j],matrix[L-1][j],labels);}
				}
			}
		}
	}

	for (int i=0; i<L; i++)
	{
		for (int j=0; j<L; j++)
		{
			if (matrix[i][j]==0) {continue;}
			matrix[i][j]=find(matrix[i][j],labels);
		}
	}
}

void HK::print_cluster()
{
	ofstream file;
	file.open("cluster.p");

	file<<"set terminal png"<<endl;
	file<<"set output 'cluster.png'"<<endl;
	file<<"set key off"<<endl;
	file<<"set xrange [0:500]"<<endl;
	file<<"set yrange [0:500]"<<endl;
	file<<"set style arrow 1 head filled size screen 0.03,15 ls 2"<<endl;

	for (int i=0; i<system.size(); i++)
	{
		for (int j=0; j<system.size(); j++)
		{
			if (matrix[i][j]==0) {continue;}
			
			double x=(i+1)*2.5;
			double y=(j+1)*2.5;

			double theta=system.angle(i*system.size()+j);
			double dx=cos(theta);
			double dy=sin(theta);

			file<<"set arrow from "<<x<<","<<y<<" to "<<x+dx<<","<<y+dy<<" as 1 lc "<<matrix[i][j]<<endl;
		}
	}
	file<<"plot NaN"<<endl;
	file.close();
}

void HK::clusters_labelled()
{
	ofstream file;
	file.open("labelled.p");

	file<<"set terminal png"<<endl;
	file<<"set output 'labelled.png'"<<endl;
	file<<"set key off"<<endl;
	file<<"set xrange [0:255]"<<endl;
	file<<"set yrange [0:255]"<<endl;
	file<<"set style arrow 1 head filled size screen 0.03,15 ls 2"<<endl;

	for (int i=0; i<system.size(); i++)
	{
		for (int j=0; j<system.size(); j++)
		{
			if (matrix[i][j]==0) {continue;}
			
			double x=(i+1)*3;
			double y=(j+1)*3;

			file<<"set label '"<<matrix[i][j]<<"' at "<<x<<","<<y<<endl;
		}
	}
	file<<"plot NaN"<<endl;
	file.close();
}

int HK::cluster_size(int label)
{
	int count=0;
	for (int i=0; i<L; i++)
	{
		for (int j=0; j<L; j++)
		{
			if (matrix[i][j]==label) {count++;}
		}
	}
	return count;
}

int HK::cluster_count()
{
	int max_label=0;
	for (int i=1; i<N_Spins; i++)
	{
		if (cluster_size(i)!=0) {max_label++;}
	}
	return max_label;
}

int HK::max_label()
{
	int max_label=0;
	for (int i=1; i<N_Spins; i++)
	{
		if (labels[i]>max_label) {max_label=labels[i];}
	}

	return max_label;
}

std::vector<double> HK::mean_distance_to_surface(int label)
{
	std::vector<int> distance;
	std::vector<xy> coord;
	std::vector<double> distr(2,0.0);

	for (int i=0; i<L; i++)
	{
		for (int j=0; j<L; j++)
		{
			if (matrix[i][j]==label)
			{
				xy temp;
				temp.x=(double) i; temp.y=(double) j;
				coord.push_back(temp);
			}
		}
	}

	double NC= (double) coord.size();

	for (int n=0; n<coord.size(); n++)
	{
		int i=coord[n].x; int j = coord[n].y;
		int nc=0.0;

		if (i== L-1 or i==0.0)
		{
			distance.push_back(0);
			continue;
		}
		else {nc += system.occupied(L*(i+1)+j)+system.occupied(L*(i-1)+j);}

		if (j==L-1 or j==0)
		{
			distance.push_back(0);
			continue;
		}
		else {nc += system.occupied(L*i+j+1)+system.occupied(L*i+j-1);}

		if (nc != 4) {distance.push_back(0); continue;}

		int D=2*L;
		for (int m=0; m<coord.size(); m++)
		{

			int s=coord[m].x; int t=coord[m].y;

			int nc2=0.0;
			if (s==L-1 or s==0) {nc2=4;}
			else {nc2 += system.occupied(L*(s+1)+t)+system.occupied(L*(s-1)+t);}

			if (t==L-1 or t==0) {nc2=4;}
			else {nc2 += system.occupied(L*s+t+1)+system.occupied(L*s+t-1);}

			int temp=abs(i-s)+abs(j-t);
			if (temp<D and nc2!=4)
			{
				D=temp;
			}

		}
		
		if (D==L) {std::cout<<"DA FUCK?"<<std::endl;}
		distance.push_back(D);
	}

	int sum=0;
	for (int u=0; u<distance.size(); u++)
	{
		sum+=distance[u];
	}
	
	double mean = (double) sum / NC;

	double sstd=0.0;
	for (int e=0; e<distance.size(); e++)
	{
		sstd+=(((double) distance[e])-mean)*(((double) distance[e])-mean);
	}
	
	double std=sqrt(sstd/(NC-1.0));
	distr[0]=mean; distr[1]=std;
	
	return distr;
}

std::vector<double> HK::mean_distance_to_surface_periodic(int label)
{
	std::vector<int> distance;
	std::vector<xy> coord;
	std::vector<double> distr(2,0.0);

	for (int i=0; i<L; i++)
	{
		for (int j=0; j<L; j++)
		{
			if (matrix[i][j]==label)
			{
				xy temp;
				temp.x=(double) i; temp.y=(double) j;
				coord.push_back(temp);
			}
		}
	}

	double NC= (double) coord.size();

	for (int n=0; n<coord.size(); n++)
	{
		int i=coord[n].x; int j = coord[n].y;
		int nc=0.0;

		if (i==0.0) {nc += system.occupied((i+1)*L+j)+system.occupied((L-1)*L+j);}
		else if (i==L-1.0) {nc += system.occupied(j)+system.occupied((i-1)*L+j);}
		else {nc += system.occupied((i+1)*L+j)+system.occupied((i-1)*L+j);}

		if (j==0.0) {nc += system.occupied(L*i+j+1)+system.occupied(i*L+L-1);}
		else if (j==L-1.0) {nc += system.occupied(i*L)+system.occupied(L*i+j-1);}
		else {nc += system.occupied(i*L+j+1)+system.occupied(L*i+j-1);}

		if (nc != 4) {distance.push_back(0); continue;}

		int D=2*L;
		for (int m=0; m<coord.size(); m++)
		{

			int s=coord[m].x; int t=coord[m].y;

			int nc2=0.0;

			if (s==0) {nc2 += system.occupied(L*(s+1)+t)+system.occupied(L*(L-1)+t);}
			else if (s==L-1) {nc2 += system.occupied(t)+system.occupied(L*(s-1)+t);}
			else {nc2 += system.occupied(L*(s+1)+t)+system.occupied(L*(s-1)+t);}

			if (t==0) {nc2 += system.occupied(L*s+t+1)+system.occupied(L*s+L-1);}
			else if (t==L-1) {nc2 += system.occupied(L*s)+system.occupied(L*s+t-1);}
			else {nc2 += system.occupied(L*s+t+1)+system.occupied(L*s+t-1);}

			int temp=abs(i-s)+abs(j-t);
			if (temp<D and nc2!=4)
			{
				D=temp;
			}

		}
		
		distance.push_back(D);
	}

	int sum=0;
	for (int u=0; u<distance.size(); u++)
	{
		sum+=distance[u];
	}
	
	double mean = (double) sum / NC;

	double sstd=0.0;
	for (int e=0; e<distance.size(); e++)
	{
		sstd+=(((double) distance[e])-mean)*(((double) distance[e])-mean);
	}
	
	double std=sqrt(sstd/(NC-1.0));
	distr[0]=mean; distr[1]=std;
	
	return distr;
}

double HK::cluster_energy(int label)
{
	double Energy=0.0;

	for (int i=0; i<system.size(); i++)
	{
		for (int j=0; j< system.size(); j++)
		{
			if (matrix[i][j]==label)
			{
				Energy+=system.H_local(i*system.size()+j);
			}
		}
	}
	return Energy;
}

double HK::cluster_energy_periodic(int label)
{
	double Energy=0.0;

	for (int i=0; i<system.size(); i++)
	{
		for (int j=0; j< system.size(); j++)
		{
			if (matrix[i][j]==label)
			{
				Energy+=system.H_local(i*L+j);
			}
		}
	}
	return Energy;
}

double HK::surface_members(int label)
{
	int on_surface=0;
	for (int i=0; i<system.size(); i++)
	{
		for (int j=0; j<system.size(); j++)
		{
			if (matrix[i][j]==label)
			{
				int neighbors;
				if (i==0)
				{
					on_surface++;
				}
				else if (j==0)
				{
					on_surface++;
				}
				else if (i==system.size())
				{
					on_surface++;
				}
				else if (j==system.size())
				{
					on_surface++;
				}
				else
				{
					neighbors=system.occupied(L*(i-1)+j)+system.occupied(L*(i+1)+j)+system.occupied(L*i+j-1)+system.occupied(L*i+j+1);
					if (neighbors<4)
					{
						on_surface++;
					}
				}
			}
		}
	}
	return on_surface;
}

double HK::circularity(int label)
{
	int L=system.size();
	double perimeter=0.0;
	double n=0.0;

	for (int i=0; i<L; i++)
	{
		for (int j=0; j< L; j++)
		{
			if (matrix[i][j]==label)
			{
				int neigh=0;
				if (i!=0) {neigh+=system.occupied(L*(i-1)+j);}
				if (i!=L-1) {neigh+=system.occupied(L*(i+1)+j);}
				if (j!=0) {neigh+=system.occupied(L*i+j-1);}
				if (j!=L-1) {neigh+=system.occupied(L*i+j+1);}

				if (neigh==4) {n+=1.0;}
				else {perimeter+=1.0;}
			}
		}
	}
	return perimeter/n;
}

std::vector<double> HK::principle_moments(int label)
{
	std::vector<double> dummy(2,0.0);
	return dummy; //This is just here to quite the compiler warning
	// std::vector<double> moments;
	// std::vector<xy> coord;

	// for (int i=0; i<L; i++)
	// {
	// 	for (int j=0; j<L; j++)
	// 	{
	// 		if (matrix[i][j]==label)
	// 		{
	// 			xy temp;
	// 			temp.x=(double) i; temp.y=(double) j;
	// 			coord.push_back(temp);
	// 		}
	// 	}
	// }

	// double NL=cluster_size(label);
	// double C=1.0/(2.0*((double) NL)*((double) NL));

	// double Sxx=0.0,
	// 	  Sxy=0.0,
	// 	  Syx=0.0,
	// 	  Syy=0.0;

	// for (int n=0; n<coord.size(); n++)
	// {
	// 	for (int m=0; m<coord.size(); m++)
	// 	{
	// 		Sxx+=C*(coord[n].x-coord[m].x)*(coord[n].x-coord[m].x);
	// 		Sxy+=C*(coord[n].x-coord[m].x)*(coord[n].y-coord[m].y);
	// 		Syx+=C*(coord[n].y-coord[m].y)*(coord[n].x-coord[m].x);
	// 		Syy+=C*(coord[n].y-coord[m].y)*(coord[n].y-coord[m].y);
	// 	}
	// }

	// double Trace = Sxx + Syy;
	// double Det = Sxx*Syy-Sxy*Syx;

	// double lambda1=0.5*Trace+sqrt(0.25*Trace*Trace - Det);
	// double lambda2=0.5*Trace-sqrt(0.25*Trace*Trace - Det);

	// if (lambda1<=lambda2)
	// {
	// 	moments.push_back(lambda1);
	// 	moments.push_back(lambda2);
	// }
	// else if(lambda1>lambda2)
	// {
	// 	moments.push_back(lambda2);
	// 	moments.push_back(lambda1);
	// }

	// return moments;
}
