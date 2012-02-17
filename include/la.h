/*
 * la.h
 */
#ifndef LA_H
#define LA_H

#include <vector>
#include "types.hpp"

namespace la {

struct Direction {
    std::vector<types::Real> incrementVector;
    types::Uint incrementCount;
};

class Cell {
    public:
        Cell() {};
        Cell(const std::vector<Direction> &directions);
        std::vector< types::Real > vector(types::Uint n) { return vectors[n]; }
    private:    
        std::vector< std::vector<types::Real> > vectors;
};


struct Grid {
    std::vector<Direction> directions;
    std::vector<types::Real> originVector;
    std::vector<types::Real> data;

    std::vector<Direction> getDirections() const { return directions;}
    std::vector<types::Real> getOriginVector() const { return originVector; }
    const std::vector<types::Real> & getData() const { return data; }
    Cell cell() const;

    const Grid & operator+=(const Grid &g);
    types::Uint countPoints() const;
    void sumXY(std::vector<types::Real>& reduced) const;
    void averageXY(std::vector<types::Real>& reduced) const;
    void squareValues();
    bool getNearestIndices(std::vector<types::Real>& coordinates, std::vector<types::Uint>& indices) const;
    types::Real getNearestDataPoint(std::vector<types::Real>& coordinates) const;
    types::Real getNearestDataPoint(types::Real x, types::Real y, types::Real z) const;
    types::Real getDataPoint(types::Uint x, types::Uint y, types::Uint z) const;
    types::Real getDataPoint(const std::vector<types::Uint>& indices) const; 
    bool checkRange(const std::vector<types::Uint>& indices) const;
    bool checkDimension(types::Uint size) const;
    bool hasSameGrid(const Grid &g) const;
    void resize(const std::vector<types::Uint>& incrementCounts);
    void resize(types::Uint nX, types::Uint nY, types::Uint nZ);
    private:
    void copyRecursive( 
        types::Uint directionIndex,
        std::vector<types::Real>::const_iterator &oldIt, 
        std::vector<types::Real>::iterator &newIt,  
        const std::vector<types::Uint> &newCounts);
};

}

#endif
