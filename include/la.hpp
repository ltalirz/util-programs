/*
 * la.h
 */
#ifndef LA_H
#define LA_H

#include <vector>
#include "types.hpp"

namespace formats {class Cube; }

namespace la {


class Direction {
    std::vector<types::Real> incrementVector;
    types::Uint incrementCount;

public:
    Direction() {};
    Direction(std::vector<types::Real> v, types::Uint c) :
        incrementVector(v), incrementCount(c) {};
    /**
     * Adapts Direction for given stride.
     */
    void stride(types::Uint s);
    void scaleVector(types::Real factor);
    types::Uint getIncrementCount() const { 
        return incrementCount;}
    const std::vector<types::Real> & getIncrementVector() const { 
        return incrementVector;}
    bool checkRange(types::Uint index) const;
    
    friend bool operator==(const Direction& d1, const Direction& d2);
    friend class formats::Cube;
};

class Cell {
    std::vector< std::vector<types::Real> > vectors;
public:
    Cell() {};
    Cell(const std::vector<Direction> &directions);
    std::vector< types::Real > vector(types::Uint n) {
        return vectors[n];
    }
};


struct Grid {
    std::vector<Direction> directions;
    std::vector<types::Real> originVector;
    std::vector<types::Real> data;

    std::vector<Direction> getDirections() const {
        return directions;
    }
    std::vector<types::Real> getOriginVector() const {
        return originVector;
    }
    const std::vector<types::Real> & getData() const {
        return data;
    }
    Cell cell() const;

    const Grid & operator+=(const Grid &g);
    const Grid & operator*=(types::Real x);
    void squareValues();
    void sqrt();
    void abs();
    types::Uint countPoints() const;

    void sumXY(std::vector<types::Real>& reduced) const;
    void zPlane(types::Uint index, std::vector<types::Real> &plane);
    void averageXY(std::vector<types::Real>& reduced) const;
    /**
     * In lack of a proper resampling method. No new values are calculated
     */
    void stride(types::Uint, types::Uint, types::Uint);
    void stride(std::vector<types::Uint> stride);
    void resize(const std::vector<types::Uint>& incrementCounts);
    void resize(types::Uint nX, types::Uint nY, types::Uint nZ);

    bool getNearestIndices(std::vector<types::Real>& coordinates, std::vector<types::Uint>& indices) const;
    types::Real getNearestDataPoint(std::vector<types::Real>& coordinates) const;
    types::Real getNearestDataPoint(types::Real x, types::Real y, types::Real z) const;
    types::Real getDataPoint(types::Uint x, types::Uint y, types::Uint z) const;
    types::Real getDataPoint(const std::vector<types::Uint>& indices) const;

    bool checkRange(const std::vector<types::Uint>& indices) const;
    bool checkDimension(types::Uint size) const;
    bool hasSameGrid(const Grid &g) const;

private:
    void copyRecursive(
        types::Uint directionIndex,
        std::vector<types::Real>::const_iterator &oldIt,
        std::vector<types::Real>::iterator &newIt,
        const std::vector<types::Uint> &newCounts);
    void strideRecursive(
        types::Uint directionIndex,
        std::vector<types::Real>::const_iterator &oldIt,
        std::vector<types::Real>::iterator &newIt,
        const std::vector<types::Uint> &newCounts,
        const std::vector<types::Uint> &strides,
        const std::vector<types::Uint> &rests);
};

}

#endif
