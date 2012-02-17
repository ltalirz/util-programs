/*
 * fundamental.h
 */
#ifndef ATOMISTIC_FUNDAMENTAL_H
#define ATOMISTIC_FUNDAMENTAL_H

#include "types.hpp"
#include "la.h"

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
struct EnergyLevels {
    static const types::Uint DEFAULT_N = 100;
    std::vector<types::Real> levels;
    types::Real fermi;

    EnergyLevels & operator*=(types::Real factor);
    EnergyLevels(std::vector<types::Real> levels, types::Real fermi);
    EnergyLevels() {};
    void sort();
    types::Uint count() const;
    void print() const;
    void shift(types::Real deltaE);
    void setFermiZero();
    types::Real getLevel(types::Uint i) const;
};

/**
 * A whole spectrum, containing several sets of energy levels
 * (spins, kpoints)
 */
struct Spectrum {
    static const types::Uint DEFAULT_N_SPINS = 2;
    std::vector<EnergyLevels> spins;

    Spectrum & operator*=(types::Real factor);
    bool readFromCp2k(types::String filename);
    bool readFromCp(types::String filename);
    void print() const;
    void shift(types::Real deltaE);
    void setFermiZero();
    EnergyLevels sumSpins() const;
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
    Atom(std::vector<types::Real> coordinates): coordinates(coordinates){};
    Atom(){ }
    const std::vector<types::Real> &getCoordinates () const {return coordinates;}
    types::Uint getNumber() const { return number; }
    types::Real getCharge() const { return charge; }
    types::String getSymbol() const { return symbol; }
};

}
#endif
