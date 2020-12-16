#include "xy.hpp"

class HK
{
	int L;
	int N_Spins;
	lattice system;
	std::vector<int> labels;
	std::vector<std::vector<int> > matrix;

public:

	HK(lattice init);
	void Find_Cluster();
	void Find_Cluster_periodic();
	void print_cluster();
	void clusters_labelled();

	int cluster_size(int label);
	int cluster_count();
	int max_label();

	double cluster_energy(int label);
	double cluster_energy_periodic(int label);
	double surface_members(int label);
	double circularity(int label);
	double distance_to_surface_periodic(int label);

	std::vector<double> principle_moments(int label);
	std::vector<double> mean_distance_to_surface(int label);
	std::vector<double> mean_distance_to_surface_periodic(int label);

};