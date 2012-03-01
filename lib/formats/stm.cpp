#include "formats/cube.hpp"
#include "formats/stm.hpp"
#include "types.hpp"
#include "io.hpp"
#include <cmath>

#include <boost/format.hpp>

namespace formats {
namespace stm {

using namespace types;

STS2d::STS2d(
        const std::list< formats::WfnCube > &cubes,
        types::Real height,
        types::Real eMin,
        types::Real eMax,
        types::Real deltaE,
        types::Real broadening){

    std::list< formats::WfnCube >::const_iterator
        cubeIt = cubes.begin(),
    cubeEnd = cubes.end();


    // Get atoms and proper grid
    this->fileName = cubeIt->fileName;
    this->readCubeFile();

    // Get z index of STS plane
    std::vector<types::Real> tempvector;
    tempvector.push_back(0);
    tempvector.push_back(0);
    tempvector.push_back(this->topZCoordinate() + height);
    std::vector<types::Uint> indices;
    this->grid.getNearestIndices(tempvector, indices);
    types::Uint zIndex = indices[2];
    std::cout << "STS will be performed at z-index " 
        << zIndex << "\n";

    this->title  = "STS data (z = energy)";
    this->title += str(boost::format(" at z-plane %3d") % zIndex);
    this->description  = str(boost::format(
                "Range [%1.3d V,%1.3d V], delta-e %1.4d, sigma %1.4d")
                % eMin % eMax % deltaE % broadening);


    // Adjust z dimension for energy
    Uint newIncrementCount =(unsigned int) ((eMax - eMin)/deltaE) + 1;
    // Want to keep old z extent
    this->grid.directions[2].scaleVector(
            Real(grid.directions[2].getIncrementCount())/
            Real(newIncrementCount));
    this->grid.directions[2] = 
        la::Direction(grid.directions[2].getIncrementVector(),
                newIncrementCount);
    this->grid.originVector[2] = eMin;
    grid.data = std::vector<Real>(grid.countPoints(), 0.0);
    this->broadening = broadening;

    // Get all necessary planes
    formats::WfnCube tempCube;
    while(cubeIt != cubeEnd){
        std::cout<< "Processing " << cubeIt->fileName << "\n";
        tempCube = *cubeIt;
        tempCube.readCubeFile();
        tempCube.squareValues();
        std::vector<Real> tempPlane;
        tempCube.getZPlane(zIndex, tempPlane);
        
        addLevel(tempPlane, cubeIt->energy); 
        ++cubeIt;
    }
        
    std::cout<< "Done processing cube files...\n\n";

    // Normalize sum to 1
    grid *= 1.0/grid.sum();
}



void STS2d::addLevel(const std::vector<types::Real> &plane,
               types::Real energy){
    // Get energy grid
    Real eMin = getEMin();
    Real deltaE = getDeltaE();
    Real eMax = getEMax();
    Uint nEnergies = grid.directions[2].getIncrementCount();

    // Gaussian stuff
    // \$ \frac{1}{\sigma \sqrt{2\pi} e^{-\frac{(x-\mu)^2}{2\sigma^2}} \$
    // = a e^{c(x-b)^2}
    types::Real a = 1.0/(broadening * std::sqrt(2 * M_PI));
    types::Real c = -1.0/(2 * broadening * broadening);
   
     std::vector<Real>::const_iterator 
         planeIt = plane.begin(),
         planeEnd = plane.end();
     std::vector<Real>::iterator 
         dataIt = grid.data.begin();
     while(planeIt != planeEnd){
         for(Uint z = 0; z < nEnergies; ++z){
             Real e = eMin + deltaE * z; 
             *dataIt += *planeIt * a * std::exp(c * (energy - e)*(energy - e));
             ++dataIt;
         }
         ++planeIt;
     }


}

Real STS2d::getEMin(){
    return grid.originVector[2];
}

Real STS2d::getEMax(){
    return (getEMin() +
           getDeltaE() * (grid.directions[2].getIncrementCount() - 1));
}

Real STS2d::getDeltaE(){
    return grid.directions[2].getIncrementVector()[2];
}

void StmCube::setIsoLevel(types::Real isoValue){
    getZIsoSurface(isoValue, stm);
    this->isoLevel = isoValue;
}

bool StmCube::writeIgorFile(String fileName) const {
    Uint nX = grid.directions[0].getIncrementCount();
    Uint nY = grid.directions[1].getIncrementCount();
    std::vector<Real>::const_iterator stmIt = stm.begin();

    Stream result = "";
    for(Uint x = 0; x<nX; ++x){
        for(Uint y=0; y<nY; ++y){
            result += str(boost::format("%12.6e") % *stmIt);
            result += " ";
            ++stmIt;
        }
        result += "\n";
    }

    return io::writeStream(fileName, result);
}


}
}
