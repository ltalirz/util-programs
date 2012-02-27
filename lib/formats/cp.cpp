#include "formats/cp.hpp"
#include "atomistic.hpp"
#include "types.hpp"
#include "io.hpp"

#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace formats {
namespace cp {

typedef boost::sregex_iterator RegIt;
typedef boost::smatch Match;
typedef const boost::regex Regex;

using boost::lexical_cast;
using namespace types;
namespace at = atomistic;

Spectrum & Spectrum::operator*=(Real factor){
    for(std::vector<at::EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        *it *= factor;
    }

    return *this;
}


void Spectrum::shift(Real deltaE) {
    for(std::vector<at::EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->shift(deltaE);
    }
}

void Spectrum::setFermiZero() {
    for(std::vector<at::EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->setFermiZero();
    }
}

at::EnergyLevels Spectrum::sumSpins() const {
    at::EnergyLevels e = at::EnergyLevels();
    // Spins may have different Fermi. The new Fermi is the average Fermi of all spins
    // and energy levels are shifted accordingly
    Uint counter = 0;
    Real fermiSum = 0;

    for(std::vector<at::EnergyLevels>::const_iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        e.join(*it);
        ++counter;
    }

    return e;
}

void Spectrum::print() const {
    for(std::vector<at::EnergyLevels>::const_iterator it = this->spins.begin(); it != this->spins.end(); it++) {
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

        at::EnergyLevels energyLevels = at::EnergyLevels(levels, 0);
        this->spins.push_back(energyLevels);

        ++occIt;
    }

    return true;
}


}
}
