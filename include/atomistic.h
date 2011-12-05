/*
 * classes.h
 */
#ifndef CLASSES_H
#define CLASSES_H

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

namespace atomistic {

typedef unsigned int Uint;
typedef int Int;
typedef double Real;
typedef std::string String;

struct EnergyLevels {
    static const Uint DEFAULT_N = 100;
    std::vector<Real> levels;
    Real fermi;

    EnergyLevels(std::vector<Real> levels, Real fermi);
    void sort();
    Uint count() const;
    void print() const;
};

struct Spectrum {
    static const Uint DEFAULT_N_SPINS = 2;
    std::vector<EnergyLevels> spins;

    void readFromCp2k(String filename);
    void print() const;
};

struct Atom {
    std::vector<Real> coordinates;
    Uint number;
    Real charge;
    String symbol;

    Atom(std::vector<Real> coordinates, Uint number, Real charge): coordinates(coordinates), number(number), charge(charge) {};
    Atom(std::vector<Real> coordinates): coordinates(coordinates){};
    Atom(){ }
    void  print() const;
};

struct FormatXyz {
    std::vector<Atom> atoms;
};

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
};

struct Cube {
	std::vector<Atom> atoms;
	Grid grid;
	std::vector<char> title;
	std::vector<char> description;

	void readCubeFile(String filename);
	void print() const;
	Uint countAtoms() const;

};

}
#endif
