/**
 * @file fundamental.hpp
 * Some classes describing fundamental atomistic concepts.
 */
#ifndef ATOMISTIC_FUNDAMENTAL_HPP
#define ATOMISTIC_FUNDAMENTAL_HPP

#include "types.hpp"
#include "la.hpp"

namespace atomistic {


//class EigenState {
//    types::Real energy;
//    WfnCube cube;
//
//    types::Real getEnergy() const { return energy;}
//    const WfnCube & getCube() const { return waveFunction;}
//};

/**
 * An array of energy levels.
 */
class EnergyLevels {
    std::vector<types::Real> levels;
    types::Real fermi;

public:
    EnergyLevels(std::vector<types::Real> levels, types::Real fermi);
    EnergyLevels() {};
    void join(EnergyLevels e2);
    
    
    EnergyLevels & operator*=(types::Real factor);
    void sort();
    types::Uint count() const;
    types::Uint countOccupied() const;
    void print() const;
    void shift(types::Real deltaE);
    void setFermiZero();

    const std::vector<types::Real> & getLevels() const { return levels; };
    void setLevels(const std::vector<types::Real> l) { levels = l; };
    types::Real getFermi() const { return fermi; }
    void setFermi(types::Real f) { fermi = f; }
    bool rangeCheck(types::Uint i) const;
    types::Real getLevel(types::Uint i) const {
        if(rangeCheck(i)) return levels[i-1]; 
        else return false;
    }
    void setLevel(types::Uint i, types::Real value) { 
        if(rangeCheck(i)) levels[i-1] = value; 
    }
};

/**
 * Representation of an atomic core in an atomistic simulation.
 */
struct Atom {
    std::vector<types::Real> coordinates;
    types::Uint number;
    types::Real charge;
    types::String symbol;

    Atom(std::vector<types::Real> coordinates,
         types::Uint number,
         types::Real charge) :
        coordinates(coordinates), number(number), charge(charge) {};
    Atom(std::vector<types::Real> coordinates): coordinates(coordinates) {};
    Atom() { }

    const std::vector<types::Real> &getCoordinates () const {
        return coordinates;
    }
    types::Uint getNumber() const {
        return number;
    }
    types::Real getCharge() const {
        return charge;
    }
    types::String getSymbol() const {
        return symbol;
    }

    friend class Cube;
};

}
#endif
