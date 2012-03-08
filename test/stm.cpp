#include <iostream>
#include <string>
#include <vector>

#include "formats/stm.hpp"
#include "types.hpp"

namespace stm = formats::stm;
using namespace types;
/**
 Generates isosurface and aftwards interpolates
 values at the coordinates of the isosurface
 to test if isosurface value is reproduced.

 For my example I found that the interpolated values
 differed from the isosurface value by a factor of up to 2.
 If we want to improve this, we need to change from linear
 to a better interpolation.
*/
int main(){
    // Read STM cube file, with profile
    stm::StmCube stmCube = stm::StmCube();
    stmCube.readCubeFile("bias_-1.000.cube");
    stmCube.readIgorFile("bias_-1.000.cube.1.0e-07.igor");

    // Use stscube for interpolation
    stm::StsCube stsCube = stm::StsCube();
    stsCube.setStm(stmCube);
    std::vector<Real> interpolatedPlane;
    stsCube.interpolateOnZProfile(stmCube.grid, interpolatedPlane);
   
    // Use stmCube for writing 
    stmCube.setZProfile(interpolatedPlane);
    stmCube.writeIgorFile("out.igor"); // Should be all 1.0e-7


    return 0;
}

