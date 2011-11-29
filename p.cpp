#include <iostream>
#include <iterator>
//#include <boost/program_options.hpp>

#include "atomistic.h"

namespace at = atomistic;

void readFromCp2k(){
	at::Spectrum spectrum = at::Spectrum();
	spectrum.readFromCp2k("stm.out");
	spectrum.print();
}

void readCubeFile(){
	at::Cube cube = at::Cube();
	cube.readCubeFile("2MOL-WFN_00122_1-1_23.cube");
}

int main(int ac, char* av[]){

	readCubeFile();

	 return 0;
}
