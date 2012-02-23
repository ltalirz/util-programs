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

    

Real EnergyLevels::getLevel(Uint i) const {
    if( i > 0 && i <= this->levels.size() ){
        return this->levels[i-1];
    } else {
        throw std::range_error("Level index out of bounds.");
        return 0;
    }
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


Spectrum & Spectrum::operator*=(Real factor){
    for(std::vector<EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        *it *= factor;
    }

    return *this;
}


void Spectrum::shift(Real deltaE) {
    for(std::vector<EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->shift(deltaE);
    }
}

void Spectrum::setFermiZero() {
    for(std::vector<EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->setFermiZero();
    }
}

EnergyLevels Spectrum::sumSpins() const {
    EnergyLevels e = EnergyLevels();
    // Spins may have different Fermi. The new Fermi is the average Fermi of all spins
    // and energy levels are shifted accordingly
    Uint counter = 0;
    Real fermiSum = 0;

    for(std::vector<EnergyLevels>::const_iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        e.join(*it);
        ++counter;
    }

    return e;
}

void Spectrum::print() const {
    for(std::vector<EnergyLevels>::const_iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->print();
    }
}

bool Spectrum::readFromCp(String filename) {
    String content;
    io::readFile(filename, content);

    RegIt occIt(
        content.begin(),
        content.end(),
        boost::regex("Eigenvalues (eV), kp =   1 , spin =  ([\\.\\-\\d\\s]*)"));
//    RegIt fermiIt(
//        content.begin(),
//        content.end(),
//        boost::regex("Fermi Energy \\[eV\\] :([\\.\\-\\d\\s]*)"));
//    RegIt unoccIt(
//        content.begin(),
//        content.end(),
//        boost::regex("Eigenvalues of the unoccupied subspace spin.*?iterations([\\.\\-\\d\\s]*)"));
    RegIt dummyIt;


    Regex fermiRegex("\\-?\\d\\.\\d{6}");
    Regex levelRegex("\\-?\\d+\\.\\d+");

    // Iterate over spins
    while(occIt != dummyIt) {
        String levelData = occIt->str();
        std::vector<Real> levels;
        RegIt levelIt(levelData.begin(), levelData.end(), levelRegex);
        while(levelIt != dummyIt) {
            levels.push_back(lexical_cast<Real>(levelIt->str()));
            ++levelIt;
        }

//        Match fermiMatch;
//        regex_search(fermiIt->str(), fermiMatch, fermiRegex);
//        Real fermi = lexical_cast<Real>(fermiMatch);
//        // fermi given in eV, want Ha = 2 Ry
//        fermi *= GSL_CONST_MKSA_ELECTRON_VOLT / 
//            (2* GSL_CONST_MKSA_RYDBERG);

        EnergyLevels energyLevels = EnergyLevels(levels, 0);
        this->spins.push_back(energyLevels);

        ++occIt;
    }

    return true;
}

bool Spectrum::readFromCp2k(String filename) {
    String content;
    io::readFile(filename, content);

    RegIt occIt(
        content.begin(),
        content.end(),
        boost::regex("Eigenvalues of the occupied subspace spin([\\.\\-\\d\\s]*)"));
    RegIt fermiIt(
        content.begin(),
        content.end(),
        boost::regex("Fermi Energy \\[eV\\] :([\\.\\-\\d\\s]*)"));
    RegIt unoccIt(
        content.begin(),
        content.end(),
        boost::regex("Eigenvalues of the unoccupied subspace spin.*?iterations([\\.\\-\\d\\s]*)"));
    RegIt dummyIt;


    Regex fermiRegex("\\-?\\d\\.\\d{6}");
    Regex levelRegex("\\-?\\d\\.\\d{8}");

    // Iterate over spins
    while(occIt != dummyIt) {
        String levelData = occIt->str().append(unoccIt->str());
        std::vector<Real> levels;
        RegIt levelIt(levelData.begin(), levelData.end(), levelRegex);
        while(levelIt != dummyIt) {
            levels.push_back(lexical_cast<Real>(levelIt->str()));
            ++levelIt;
        }

        Match fermiMatch;
        regex_search(fermiIt->str(), fermiMatch, fermiRegex);
        Real fermi = lexical_cast<Real>(fermiMatch);
        // fermi given in eV, want Ha = 2 Ry
        fermi *= atomistic::units::eV / atomistic::units::Ha;

        EnergyLevels energyLevels = EnergyLevels(levels, fermi);
        this->spins.push_back(energyLevels);

        ++occIt;
        ++unoccIt;
        ++fermiIt;
    }

    return true;
}

}

