#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>

#include "atomistic/fundamental.hpp"
#include "atomistic/units.hpp"
#include "io.hpp"
#include "la.hpp"
#include "types.hpp"


namespace atomistic {

typedef boost::sregex_iterator RegIt;
typedef boost::smatch Match;
typedef const boost::regex Regex;

using boost::lexical_cast;
using namespace types;


EnergyLevels::EnergyLevels(std::vector<types::Real> levels, types::Real fermi){
    this->levels = levels;
    this->fermi = fermi;
    this->sort();
}


EnergyLevels & EnergyLevels::operator*=(Real factor){
    for(std::vector<Real>::iterator it = this->levels.begin(); it != this->levels.end(); it++) {
        *it *= factor;
    }
    fermi *= factor;
    return *this;
}



void EnergyLevels::sort() {
    std::sort(this->levels.begin(), this->levels.end());
}

Uint EnergyLevels::count() const {
    return this->levels.size();
}

Uint EnergyLevels::countOccupied() const {
    Uint n = 0;
    for(std::vector<Real>::const_iterator it = this->levels.begin();
            it != this->levels.end(); it++) {
        if(*it < fermi) ++n;
    }
    return n;
}


void EnergyLevels::join(EnergyLevels e2){
        e2.shift(fermi - e2.fermi);
        levels.insert( levels.end(), e2.levels.begin(), e2.levels.end());
        this->sort();
}


bool EnergyLevels::rangeCheck(Uint i) const {
    if( i <= this->levels.size()) return true;
#ifdef ATOMISTIC_FUNDAMENTAL_STRICT
    else throw std::range_error("Level index out of bounds.");
#else
    else return false;
#endif
}

void EnergyLevels::print() const {
    std::cout << "Fermi energy: " << this->fermi << std::endl;
    for(std::vector<Real>::const_iterator it = this->levels.begin(); it != this->levels.end(); it++) {
        std::cout << *it << std::endl;
    }
}

void EnergyLevels::shift(Real deltaE) {
    for(std::vector<Real>::iterator it = this->levels.begin(); it != this->levels.end(); it++) {
        *it += deltaE;
    }
    fermi += deltaE;
}

void EnergyLevels::setFermiZero() {
    this->shift(-fermi);
}



}

