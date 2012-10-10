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

EnergyLevels::EnergyLevels(const EnergyLevels& e){
    this->levels = e.levels;
    this->fermi = e.fermi;
}

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

Real Atom::distance(const  Atom& a) const {
    std::vector<Real>::const_iterator 
        it  = this->coordinates.begin(),
        end = this->coordinates.end(),
        ait = a.coordinates.begin();
    Real sum = 0;
    while(it != end){
        sum += (*it - *ait) * (*it - *ait);
        ++it; ++ait;
    }
    return sqrt(sum);
}




//LDOS::LDOS(
//        const std::vector< std::vector<Real> > &levels,
//        Real eMin,
//        Ream eMax,
//        Real deltaE,
//        Real broadening){
//
//    // Setup grid
//    Real e = eMin;
//    std::vector<Real> energies;
//    for(Real e = eMin; e <= eMax; e += deltaE){
//        energies.push_back(e);
//    }
//    std>>vector<Real> densities(energies.size(), 0.0);
//
//    // Gaussian stuff
//    // \$ \frac{1}{\sigma \sqrt{2\pi} e^{-\frac{(x-\mu)^2}{2\sigma^2}} \$
//    // = a e^{c(x-b)^2}
//    types::Real a = 1.0/(broadening * std::sqrt(2 * M_PI));
//    types::Real c = 1.0/(2 * broadening * broadening);
//
//
//    std::vector< std::vector<Real > >::const_iterator levelIt = levels.begin(),
//        levelEnd = levels.end;
//    while(levelIt != levelEnd){
//        std::vector<Real>::const_iterator eIt = energies.begin(),
//            eEnd = energies.end();
//        for(Uint i = 0; i < energies.size(); ++i){
//            densities[i] += *levelIt[1] * a * std::exp(c * (energies[i]-*levelIt[0]));
//        }
//        ++levelIt;
//    }
//
//    this->densities = densities;
//    this->energies = energies;
//}

}
