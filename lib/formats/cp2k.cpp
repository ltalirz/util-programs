#include "formats/cp2k.hpp"
#include "atomistic.hpp"
#include "types.hpp"
#include "io.hpp"

#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace formats {
namespace cp2k {

typedef boost::sregex_iterator RegIt;
typedef boost::smatch Match;
typedef const boost::regex Regex;

using boost::lexical_cast;
using namespace types;

Spectrum::Spectrum(const Spectrum& s){
    this->spins = s.spins;
}

Spectrum & Spectrum::operator*=(Real factor){
    for(std::vector<atomistic::EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        *it *= factor;
    }

    return *this;
}

void Spectrum::sort(){
    for(std::vector<atomistic::EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->sort();
    }

}

void Spectrum::shift(Real deltaE) {
    for(std::vector<atomistic::EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->shift(deltaE);
    }
}

void Spectrum::setFermiZero() {
    for(std::vector<atomistic::EnergyLevels>::iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->setFermiZero();
    }
}

atomistic::EnergyLevels Spectrum::sumSpins() const {
    atomistic::EnergyLevels e = atomistic::EnergyLevels();
    // Spins may have different Fermi. The new Fermi is the average Fermi of all spins
    // and energy levels are shifted accordingly
    Uint counter = 0;
    Real fermiSum = 0;

    for(std::vector<atomistic::EnergyLevels>::const_iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        e.join(*it);
        ++counter;
    }

    return e;
}

void Spectrum::print() const {
    for(std::vector<atomistic::EnergyLevels>::const_iterator it = this->spins.begin(); it != this->spins.end(); it++) {
        it->print();
    }
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

    bool hasUnocc = (unoccIt != dummyIt);


    Regex fermiRegex("\\-?\\d\\.\\d{6}");
    Regex levelRegex("\\-?\\d\\.\\d{8}");

    // Iterate over spins
    while(occIt != dummyIt) {
        String levelData = occIt->str();
        // Unoccupied levels are optional
        if(hasUnocc) levelData.append(unoccIt->str());
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

        atomistic::EnergyLevels energyLevels = atomistic::EnergyLevels(levels, fermi);
        this->spins.push_back(energyLevels);

        ++occIt;
        if(hasUnocc) ++unoccIt;
        ++fermiIt;
    }

    if(this->spins.size() == 0) 
        throw types::parseError() << types::errinfo_parse("No energy levels found in file.");

    return true;
}

Real Spectrum::getLevel(Uint nSpin, Uint nLevel) const {
    return spins[nSpin -1].getLevel(nLevel);
}

}
}
