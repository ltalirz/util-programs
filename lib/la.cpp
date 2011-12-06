#include "la.h"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <boost/format.hpp>
#include <eigen3/Eigen/Dense>

namespace la {

Int round(Real r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}


void Grid::printHeader() const {
    using boost::format;
    format formatter("%5i %12.6f %12.6f %12.6f \n");

    std::vector<Direction>::const_iterator it;
    for(it = directions.begin(); it!= directions.end(); ++it) {
        std::cout << formatter % it->incrementCount % it->incrementVector[0] % it->incrementVector[1] % it->incrementVector[2];
    }

}

void Grid::printData() const {
    using boost::format;
    format formatter("%13.5e");

    std::vector<Real>::const_iterator it = data.begin();
    // Fastest direction is z, stored in directions[2]
    Uint i = 1, nZ = directions[2].incrementCount;
    while(it!= data.end()) {
        std::cout << formatter % *it;
        if(i % 6 == 0) std::cout << std::endl;
        else if(i % nZ == 0) {
            std::cout << std::endl;
            i=0;
        }
        ++it;
        ++i;
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


// 3d wrapper for nd getData
Real Grid::getData(Uint x, Uint y, Uint z) const {
    std::vector<Uint> indices;
    indices.push_back(x);
    indices.push_back(y);
    indices.push_back(z);

    return this->getData(indices);
}


// This works for n dimensions,
// assuming that they are ordered by increasing fastness
Real Grid::getData(const std::vector<Uint>& indices) const {
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

// 3d wrapper for nd getData
Real Grid::getNearestData(Real x, Real y, Real z) const {
    std::vector<Real> indices;
    indices.push_back(x);
    indices.push_back(y);
    indices.push_back(z);

    return this->getNearestData(indices);
}

// This works for nd grids
Real Grid::getNearestData(std::vector<Real>& cartesianCoordinates) const {

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
    std::vector<Uint> gridStlRounded;
    
    for(std::vector<Real>::iterator it = gridStl.begin(); it != gridStl.end(); ++it) {
        gridStlRounded.push_back(round(*it));
    }

    return this->getData(gridStlRounded);
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

