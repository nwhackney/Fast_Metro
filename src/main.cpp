#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <limits>

#include "../include/gsl/gsl_rng.h"
//#include "../include/lattice.hpp"
#include "../include/simulation.hpp"
#include "../include/tinytoml-master/include/toml/toml.h"

using namespace std;

int main()
{
	simulation A;
	A.read_config();
	A.run();

	return 0;
}