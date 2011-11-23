/*
 * classes.h
 */
#ifndef CLASSES_H
#define CLASSES_H

typedef uint unsigned int;
typedef sint int;
typedef real double;
typedef text std::string;

struct EnergyLevels {
	std::vector<real> occupied;
	std::vector<real> unOccupied;
	real fermi;
}

struct Atom {
	real x,y,z;
	text symbol;
}
struct FormatXyz {
	std::vector<Atom> atoms;
}

struct Cell {
	std::vector<real>(3) a,b,c;
}

#endif
