#include <vector>
#include <math.h>
#include <iostream>


class simulation
{
	int L,
	    N,
	    Time,
	    end_time,
	    restart_t,
	    seed,
	    Num_Rep;

	double slp,
		  w,
		  J,
		  K,
		  f,
		  Tf;
	
	std::string restart,
		       out_file,
		       in_file,
	            run_type;

	gsl_rng * r;

public:
	void read_config();
	void simulated_annealing();
	void parallel_tempering();
	void run();
};