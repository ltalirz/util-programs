/*
 * classes.h
 */
#ifndef CLASSES_H
#define CLASSES_H

#include "types.h"
#include "la.h"

namespace atomistic {

/**
 * An array of energy levels.
 */
struct EnergyLevels {
    static const types::Uint DEFAULT_N = 100;
    std::vector<types::Real> levels;
    types::Real fermi;

    EnergyLevels(std::vector<types::Real> levels, types::Real fermi);
    void sort();
    types::Uint count() const;
    void print() const;
    void shift(types::Real deltaE);
    types::Real getLevel(types::Uint i) const;
};

/**
 * A whole spectrum, containing several sets of energy levels
 * (spins, kpoints)
 */
struct Spectrum {
    static const types::Uint DEFAULT_N_SPINS = 2;
    std::vector<EnergyLevels> spins;

    bool readFromCp2k(types::String filename);
    void print() const;
    void shift(types::Real deltaE);
};

/** 
 * Representation of an atomic core in an atomistic simulation.
 */
struct Atom {
    std::vector<types::Real> coordinates;
    types::Uint number;
    types::Real charge;
    types::String symbol;

    Atom(std::vector<types::Real> coordinates, types::Uint number, types::Real charge): coordinates(coordinates), number(number), charge(charge) {};
    Atom(std::vector<types::Real> coordinates): coordinates(coordinates){};
    Atom(){ }
    const std::vector<types::Real> &getCoordinates () const { return coordinates;}
    types::Uint getNumber() const { return number; }
    types::Real getCharge() const { return charge; }
    types::String getSymbol() const { return symbol; }
};

/**
 * Stores .xyz file
 */
struct FormatXyz {
    std::vector<Atom> atoms;
};

/**
 * Stores .cube file
 */
struct Cube {
	std::vector<Atom> atoms;
	la::Grid grid;
    types::Binary title;
    types::Binary description;

	bool readCubeFile(types::String filename);
	bool writeCubeFile(types::String filename) const;
	bool writeZProfile(types::String filename, types::String header) const;
	bool writeZProfile(types::String filename) const;
        void addZProfile(types::Stream &stream, types::String header) const;
	void print() const;
	void addHeader(types::Stream &stream) const;
	void addData(types::Stream &stream) const;
	types::Uint countAtoms() const;
	types::Uint countPoints() const;
	types::Real getEnergy() const;


};

/**
 * Cube file of a wave function.
 * CP2K will write the level and the spin in the title
 */
struct WfnCube : Cube {
    types::Uint spin;
    types::Uint wfn;
    types::Real energy;
    void readCubeFile(types::String filename);
};

}
#endif
