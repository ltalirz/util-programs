/*
 * la.h
 */
#ifndef LA_H
#define LA_H

#include <vector>
#include "types.h"

namespace la {

using namespace types;

struct Direction {
    std::vector<Real> incrementVector;
    Uint incrementCount;
};

struct Grid {
    std::vector<Direction> directions;
    std::vector<Real> originVector;
    std::vector<Real> data;

    void printHeader() const;
    void printData() const;
    Uint countPoints() const;
    void sumXY(std::vector<Real>& reduced) const;
    void squareValues();
    Real getNearestData(std::vector<Real>& coordinates) const;
    Real getNearestData(Real x, Real y, Real z) const;
    Real getData(Uint x, Uint y, Uint z) const;
    Real getData(const std::vector<Uint>& indices) const; 
    bool checkRange(const std::vector<Uint>& indices) const;
    bool checkDimension(Uint size) const;
};

}

#endif
