#include "la.hpp"
#include "types.hpp"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cmath>
#include <boost/format.hpp>

#include <blitz/array.h>

#include <eigen3/Eigen/Dense>


namespace la {

using namespace types;

Int round(Real r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

Cell::Cell(const std::vector<Direction> &directions) {
    for(std::vector<Direction>::const_iterator dIt = directions.begin();
            dIt != directions.end(); ++dIt) {
        Direction temp = *dIt;
        temp.scaleVector(temp.getIncrementCount());
        vectors.push_back(temp.getIncrementVector());
    }
}

Cell Grid::cell() const {
    return Cell(directions);
}

bool Grid::hasSameGrid(const Grid &g) const {
    if(this->directions == g.directions &&
            this->originVector == g.originVector)
        return true;
    else return false;
}


bool operator==(const Direction& d1, const Direction& d2) {
    return (d1.incrementVector == d2.incrementVector &&
            d1.incrementCount == d2.incrementCount);
}

bool Direction::checkRange(Uint index) const {
    if(index < incrementCount) return true;
#ifdef LA_STRICT
    else throw std::range_error("Index out of range.");
#else
    else return false;
#endif
}

void Direction::stride(types::Uint s) {
    this->scaleVector(s);
    if(incrementCount % s == 0) incrementCount = incrementCount/s;
    else incrementCount = incrementCount/s + 1;
}

void Direction::scaleVector(types::Real factor) {
    for(std::vector<types::Real>::iterator vIt = incrementVector.begin();
            vIt != incrementVector.end(); ++vIt) {
        *vIt *= factor;
    }
}

const Grid & Grid::operator+=(const Grid &g) {
    if(!this->hasSameGrid(g))
        throw types::runtimeError() <<
                                    types::errinfo_runtime("Cannot add data on different grids.");

    else if(this->data.size() == g.data.size() && this->data.size() > 0) {
        std::vector<Real>::iterator dataIt = data.begin(), dataEnd = data.end();
        std::vector<Real>::const_iterator gIt = g.data.begin();

        while(dataIt != dataEnd) {
            *dataIt += *gIt;
            ++dataIt;
            ++gIt;
        }
    }
    return *this;
}

void Grid::stride(std::vector<Uint> strides) {
    Uint dimension = strides.size();
    checkDimension(dimension);

    // Get size for preallocation
    Uint size = 1;
    std::vector<Direction>::iterator dirIt = directions.begin();
    std::vector<Uint>::iterator sIt = strides.begin();
    std::vector<Uint> rests;
    std::vector<Uint> newCounts;
    std::vector<Direction> newDirections;
    while(dirIt != directions.end()) {
        Direction d = *dirIt;
        d.stride(*sIt);
        size *= d.getIncrementCount();
        newDirections.push_back(d);
        newCounts.push_back(d.getIncrementCount());

        // Handle "rest"
        Uint r = dirIt->getIncrementCount() % *sIt;
        if(r > 0) rests.push_back(r);
        else rests.push_back(*sIt);

        ++dirIt;
        ++sIt;
    }
    std::vector<Real> newData;
    newData.resize(size);

    // Starting with direction 0, the slowest direction
    std::vector<Real>::const_iterator dataIt = data.begin();
    std::vector<Real>::iterator newIt = newData.begin();
    this->strideRecursive(0, dataIt, newIt, newCounts, strides, rests);

    // Update data and incrementCounts
    this->data = newData;
    this->directions = newDirections;
}

void Grid::strideRecursive(
    Uint directionIndex,
    std::vector<Real>::const_iterator &oldIt,
    std::vector<Real>::iterator &newIt,
    const std::vector<Uint> &newCounts,
    const std::vector<Uint> &strides,
    const std::vector<Uint> &rests) {

    Uint dimension = strides.size();

    // If directionIndex is out of bounds, we are finished
    if( directionIndex >= dimension) return;

    Uint oldCount = directions[directionIndex].getIncrementCount(),
         newCount = newCounts[directionIndex];

    Uint weight = 1;
    for(Uint i = dimension -1; i > directionIndex; --i) {
        weight *= directions[i].getIncrementCount();
    }
    Uint stride = weight * strides[directionIndex];
    Uint rest = weight * rests[directionIndex];

    // If we are at the last dimension, we copy and move the iterators
    if( directionIndex == dimension - 1) {
        for(Uint i = 0; i < newCount; ++i) {
            *newIt = *oldIt;
            ++newIt;
            oldIt += stride;
        }
    }
    else {
        for(Uint i = 0; i < newCount; ++i) {
            this->strideRecursive(directionIndex + 1, oldIt, newIt, newCounts, strides, rests);
            // The above already moved oldIt by 1xweigth
            oldIt += stride - weight;
        }
    }

    // The above is sufficient, if stride divides incrementCount
    // In this case, rest = stride
    oldIt += rest;
    oldIt -= stride;
}






/** Can cut down as well as enlarge (fill with zeros)
 *
 * Nd implementation sadly needs recursion (depth = number of dimensions).
 * I should test the speed and maybe think about a specialization for 3d.
 * Maybe one can do the recursion during compile time via templates.
 */
void Grid::resize(const std::vector<Uint>& incrementCounts) {
    Uint dimension = incrementCounts.size();
    checkDimension(dimension);

    // Get size for preallocation
    Uint size = 1;
    std::vector<Uint>::const_iterator incIt = incrementCounts.begin(),
                                      incEnd = incrementCounts.end();
    while(incIt != incEnd) {
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
    for(incIt = incrementCounts.begin(); incIt != incEnd; ++incIt, ++dirIt) {
        Direction temp = Direction(dirIt->getIncrementVector(), *incIt);
        *dirIt = temp;
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
    const std::vector<Uint> &newCounts) {

    Uint dimension = newCounts.size();

    // If directionIndex is out of bounds, we are finished
    if( directionIndex >= dimension) return;

    Uint oldCount = directions[directionIndex].getIncrementCount(),
         newCount = newCounts[directionIndex];
    Uint copyCount = std::min(oldCount, newCount),
         iterCount = std::max(oldCount, newCount);

    // If we are at the last dimension, we copy and move the iterators
    if( directionIndex == dimension - 1) {
        copy(oldIt, oldIt + copyCount, newIt);
        oldIt += oldCount;
        newIt += newCount;
    }
    // If we are not at the last dimension, we call recursively
    else {

        for(Uint i = 0; i < iterCount; ++i) {
            // If we need to copy...
            if(i < copyCount) {
                this->copyRecursive(directionIndex + 1, oldIt, newIt, newCounts);
            }
            // Else countNew or countOld is > countCopy and we just need to move iterators.
            else {
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

void Grid::abs() {
    std::vector<Real>::iterator it = data.begin(), end = data.end();
    while(it != end) {
        *it = std::abs(*it);
        ++it;
    }
}
void Grid::sqrt() {
    std::vector<Real>::iterator it = data.begin(), end = data.end();
    while(it != end) {
        *it = std::sqrt(std::abs(*it));
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
        if( !dirIt->checkRange(*indexIt) ) return false;
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
        weight *= dirIt->getIncrementCount();
        dataIt += (*indexIt) * weight;
    }

    return *dataIt;
}

// This works for nd grids
bool Grid::getNearestIndices(std::vector<Real>& cartesianCoordinates, std::vector<Uint>& indices) const {
    Uint dim = cartesianCoordinates.size();
    if (!checkDimension(dim)) return false;

    std::vector<Real> basisVectors;
    std::back_insert_iterator< std::vector<Real> > it(basisVectors);
    std::vector<Direction>::const_iterator dirIt = directions.begin(), dirEnd = directions.end();
    while(dirIt != dirEnd) {
        std::vector<Real> temp = dirIt->getIncrementVector();
        copy(temp.begin(), temp.end(), it);
        ++dirIt;
    }


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

const Grid & Grid::operator*=(types::Real x) {
    std::vector<Real>::iterator dataIt = data.begin(),
                                dataEnd = data.end();
    while(dataIt != dataEnd) {
        *dataIt *= x;
        ++dataIt;
    }
    return *this;
}

Uint Grid::countPoints() const {
    std::vector<Direction>::const_iterator it;
    Uint points = 1;
    for(it = directions.begin(); it!= directions.end(); ++it) {
        points *= it->getIncrementCount();
    }
    return points;
}




}

