#include <iostream>

#include "atomistic.h"

namespace at = atomistic;

void readFromCp2k(){
	at::Spectrum spectrum = at::Spectrum();
	spectrum.readFromCp2k("stm.out");
	spectrum.print();
}

void readCubeFile(){
	at::WfnCube cube = at::WfnCube();
	cube.readCubeFile("2MOL-WFN_00122_1-1_23.cube");
	
	std::cout << cube.grid.getData(0,0,0) << std::endl;
	std::cout << cube.grid.getData(0,1,0) << std::endl;

	std::cout << cube.grid.getNearestData(0,0,0) << std::endl;


}

int main(int ac, char* av[]){

	readCubeFile();

	 return 0;
}
