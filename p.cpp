#include <iostream>
#include <iterator>
//#include <boost/program_options.hpp>

#include "include/classes.h"

namespace at = atomistic;

int main(int ac, char* av[]){

	at::EnergyLevels levels = at::EnergyLevels();
	levels.readFromCp2k("bla");

	 return 0;
}
