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

	typedef unsigned int uint;
	typedef int sint;
	typedef double real;
	typedef std::string text;

	

	struct EnergyLevels {
		static const uint DEFAULT_N = 100;
		std::vector<real> levels;
		real fermi;
		
		void sort(){ std::sort(levels.begin(), levels.end()); }
		void readFromCp2k(std::string filename);
	};

	struct Atom {
		real x,y,z;
		text symbol;
	};

	struct FormatXyz {
		std::vector<Atom> atoms;
	};

	struct Cell {
		std::vector<real> a,b,c;
	};

}

#endif
