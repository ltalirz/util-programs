#include "la.h"
#include "types.h"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <boost/format.hpp>

#include <eigen3/Eigen/Dense>


namespace la {

using namespace types;

Int round(Real r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

Cell::Cell(const std::vector<Direction> &directions){
    for(std::vector<Direction>::const_iterator dIt = directions.begin();
            dIt != directions.end(); ++dIt){
        std::vector<types::Real> cellVector;
        for(std::vector<types::Real>::const_iterator vIt = dIt->incrementVector.begin();
                vIt != dIt->incrementVector.end(); ++vIt){
            cellVector.push_back(*vIt * dIt->incrementCount);
        }
        vectors.push_back(cellVector);
    }
}

Cell Grid::cell() const {
    return Cell(directions);
}
/**
 * 3d wrapper for general resize function
 */
void Grid::resize(types::Uint nX, types::Uint nY, types::Uint nZ){
    std::vector<Uint> tempCounts;
    tempCounts.push_back(nX);
    tempCounts.push_back(nY);
    tempCounts.push_back(nZ);
    resize(tempCounts);
}

/** Can cut down as well as enlarge (fill with zeros)
 *
 * Nd implementation sadly needs recursion (depth = number of dimensions).
 * I should test the speed and maybe think about a specialization for 3d.
 * Maybe one can do the recursion during compile time via templates.
 */
void Grid::resize(const std::vector<Uint>& incrementCounts){
    Uint dimension = incrementCounts.size();
    checkDimension(dimension);
  
    // Get size for preallocation
    Uint size = 1;
    std::vector<Uint>::const_iterator incIt = incrementCounts.begin(),
        incEnd = incrementCounts.end();
    while(incIt != incEnd){
        size *= *incIt;
	    incIt++;
    }
    std::vector<Real> newData;
    newData.resize(size);
   
    // Starting with direction 0, the slowest direction
    std::vector<Real>::const_iterator dataIt = data.begin();
    std::vector<Real>::iterator newIt = newData.begin();
    this->copyRecursive(0, dataIt, newIt, incrementCounts);
   
   // Update data and incrementCounts
   this->data = newData;
   std::vector<Direction>::iterator dirIt = directions.begin();
   for(incIt = incrementCounts.begin(); incIt != incEnd; ++incIt, ++dirIt){
      dirIt->incrementCount = *incIt;
   }
}

/** Used by resize (recursively).
 * 
 * The *last* dimension is the fast one (as in 3d cube file format)
 */
void Grid::copyRecursive( 
        Uint directionIndex,
        std::vector<Real>::const_iterator &oldIt, 
        std::vector<Real>::iterator &newIt,  
        const std::vector<Uint> &newCounts){

    Uint dimension = newCounts.size();

    // If directionIndex is out of bounds, we are finished
    if( directionIndex >= dimension) return;

    Uint oldCount = directions[directionIndex].incrementCount,
         newCount = newCounts[directionIndex];
    Uint copyCount = std::min(oldCount, newCount),
         iterCount = std::max(oldCount, newCount);
    
    // If we are at the last dimension, we copy and move the iterators
    if( directionIndex == dimension - 1){
        copy(oldIt, oldIt + copyCount, newIt);
        oldIt += oldCount;
        newIt += newCount; 
    }
    // If we are not at the last dimension, we call recursively  
    else{

        for(Uint i = 0; i < iterCount; ++i){
            // If we need to copy...
            if(i < copyCount){
                this->copyRecursive(directionIndex + 1, oldIt, newIt, newCounts);
            }
            // Else countNew or countOld is > countCopy and we just need to move iterators.
            else{
                if(newCount > i)
                    newIt += newCount;
                else
                    oldIt += oldCount;
            }
        }

    }
}

void Grid::squareValues() {
    std::vector<Real>::iterator it = data.begin(), end = data.end();
    while(it != end) {
        (*it) *= (*it);
        ++it;
    }
}

bool Grid::checkDimension(Uint size) const {
    if (size != directions.size()) {
        throw std::range_error("Number of given indices does not equal number of grid directions.");
    }

    return true;
}

bool Grid::checkRange(const std::vector<Uint>& indices) const {

    checkDimension(indices.size());
    std::vector<Uint>::const_iterator indexIt = indices.end(), indexEnd = indices.end();
    std::vector<Direction>::const_iterator dirIt = directions.begin();
    while(indexIt != indexEnd) {
	// index may range from 0..incrementCount
        if(*indexIt >= dirIt->incrementCount) {
            throw std::range_error("Given index is out of range.");
        }
        ++indexIt;
        ++dirIt;
    }

    return true;
}

Real Grid::getNearestDataPoint(std::vector<Real>& coordinates) const {
	std::vector<Uint> indices;
this->getNearestIndices(coordinates, indices);
    return this->getDataPoint(indices);
}


// 3d wrapper for nd getNearestDataPoint
Real Grid::getNearestDataPoint(Real x, Real y, Real z) const {
    std::vector<Real> coordinates;
    coordinates.push_back(x);
    coordinates.push_back(y);
    coordinates.push_back(z);

    return this->getNearestDataPoint(coordinates);
}


// 3d wrapper for nd getDataPoint
Real Grid::getDataPoint(Uint x, Uint y, Uint z) const {
    std::vector<Uint> indices;
    indices.push_back(x);
    indices.push_back(y);
    indices.push_back(z);

    return this->getDataPoint(indices);
}


// This works for n dimensions,
// assuming that they are ordered by increasing fastness
Real Grid::getDataPoint(const std::vector<Uint>& indices) const {
    checkRange(indices);
    
    // Note: end() points *after* the last item
    std::vector<Uint>::const_iterator indexIt = indices.end(), indexBegin = indices.begin();
    std::vector<Direction>::const_iterator dirIt = directions.end(), dirFirst = directions.begin();

    --indexIt;
    std::vector<Real>::const_iterator dataIt = data.begin();
    dataIt += *indexIt;
    Uint weight = 1;
    while(dirIt != dirFirst) {
        --dirIt;
        --indexIt;
        weight *= dirIt->incrementCount;
        dataIt += (*indexIt) * weight;
    }

    return *dataIt;
}

// This works for nd grids
bool Grid::getNearestIndices(std::vector<Real>& cartesianCoordinates, std::vector<Uint>& indices) const {

    std::vector<Real> basisVectors;
    std::back_insert_iterator< std::vector<Real> > it(basisVectors);
    std::vector<Direction>::const_iterator dirIt = directions.begin(), dirEnd = directions.end();
    while(dirIt != dirEnd) {
        copy(dirIt->incrementVector.begin(), dirIt->incrementVector.end(), it);
        ++dirIt;
    }

    Uint dim = cartesianCoordinates.size();
    checkDimension(dim);

    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::Map;
    Map<MatrixXd> basisMatrix(&basisVectors[0], dim, dim);
    Map<VectorXd> cartesianEig(&cartesianCoordinates[0], dim);
    VectorXd gridEig = basisMatrix.inverse() * cartesianEig;

    std::vector<Real> gridStl(gridEig.data(), gridEig.data() + gridEig.size());
    indices.clear(); 
    for(std::vector<Real>::iterator it = gridStl.begin(); it != gridStl.end(); ++it) {
        indices.push_back(round(*it));
    }

    return true;
}

void Grid::sumXY(std::vector<Real>& reduced) const {
    /** The Blitz++ way
    using namespace blitz;
    namespace t = tensor;

    Array<Real,3> dataArray(&data[0], shape(directions[0].incrementCount, directions[1].incrementCount, directions[2].incrementCount));
    // reduce second dimension
    Array<Real,2> reducedY(sum(dataArray(t::i,t::k,t::j), t::k));
    //std::cout << reduceY;
    Array<Real,1> reducedXY(sum(reducedY(t::j,t::i), t::j));
    std::cout << reduceXY;

     **/

    reduced = std::vector<Real> (directions[2].incrementCount, 0.0);
    std::vector<Real>::const_iterator itData=data.begin(), endData=data.end();
    std::vector<Real>::iterator itReduced=reduced.begin(), endReduced=reduced.end();
    // z is the fast index of the cube file, so we just need to sum
    // all of the z-compartments together
    while(itData != endData) {
        if(itReduced == endReduced) {
            itReduced = reduced.begin();
        }
        *itReduced += *itData;
        ++itData;
        ++itReduced;
    }


}

Uint Grid::countPoints() const {
    std::vector<Direction>::const_iterator it;
    Uint points = 1;
    for(it = directions.begin(); it!= directions.end(); ++it) {
        points *= it->incrementCount;
    }
    return points;
}




}

