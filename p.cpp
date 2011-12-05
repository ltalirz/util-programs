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
	cube.grid.squareValues();
	std::vector<double> line;
	cube.grid.sumXY(line);
	std::vector<double>::iterator it;
	for(it = line.begin();it!=line.end();++it){
	    std::cout << *it << std::endl;
	}
}

int main(int ac, char* av[]){

	readCubeFile();

	 return 0;
}
