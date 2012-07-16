#include "formats/espresso.hpp"
#include "atomistic.hpp"
#include "types.hpp"
#include "io.hpp"

#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace formats {
namespace espresso {

typedef boost::sregex_iterator RegIt;
typedef boost::smatch Match;
typedef const boost::regex Regex;

using boost::lexical_cast;
using namespace types;

Spectrum & Spectrum::operator*=(Real factor){
    for(std::vector<atomistic::EnergyLevels>::iterator it = this->kpoints.begin(); it != this->kpoints.end(); it++) {
        *it *= factor;
    }

    return *this;
}


void Spectrum::shift(Real deltaE) {
    for(std::vector<atomistic::EnergyLevels>::iterator it = this->kpoints.begin(); it != this->kpoints.end(); it++) {
        it->shift(deltaE);
    }
}

void Spectrum::setFermiZero() {
    for(std::vector<atomistic::EnergyLevels>::iterator it = this->kpoints.begin(); it != this->kpoints.end(); it++) {
        it->setFermiZero();
    }
}

atomistic::EnergyLevels Spectrum::mergeKpoints() const {
    atomistic::EnergyLevels e = atomistic::EnergyLevels();

    for(std::vector<atomistic::EnergyLevels>::const_iterator it = this->kpoints.begin(); it != this->kpoints.end(); it++) {
        e.join(*it);
    }

    return e;
}

void Spectrum::print() const {
    for(std::vector<atomistic::EnergyLevels>::const_iterator it = this->kpoints.begin(); it != this->kpoints.end(); it++) {
        it->print();
    }
}

bool Spectrum::readFromOutput(String filename) {
    String content;
    io::readFile(filename, content);

    RegIt kpIt(
        content.begin(),
        content.end(),
        boost::regex("k =.*? \\(ev\\):([\\.\\-\\d\\s]*)"));
    RegIt fermiIt(
        content.begin(),
        content.end(),
        boost::regex("the Fermi energy is ([\\.\\-\\d\\s]*)"));
    RegIt dummyIt;


    Regex levelRegex("\\-?\\d*\\.\\d{4}");

    Match fermiMatch;
    regex_search(fermiIt->str(), fermiMatch, levelRegex);
    Real fermi = lexical_cast<Real>(fermiMatch);
    
    // Iterate over kpoints
    while(kpIt != dummyIt) {
        String levelData = kpIt->str();
        std::cout << levelData;
        
        std::vector<Real> levels;
        RegIt levelIt(levelData.begin(), levelData.end(), levelRegex);
        
        // The first three matched numbers are the kpoint coordinates
        std::vector<Real> kVector;
        kVector.push_back(lexical_cast<Real>(levelIt->str())); ++levelIt;
        kVector.push_back(lexical_cast<Real>(levelIt->str())); ++levelIt;
        kVector.push_back(lexical_cast<Real>(levelIt->str())); ++levelIt;

        while(levelIt != dummyIt) {
            levels.push_back(lexical_cast<Real>(levelIt->str()));
            ++levelIt;
        }

        atomistic::EnergyLevels energyLevels = atomistic::EnergyLevels(levels, fermi);
        this->kpoints.push_back(energyLevels);
        ++kpIt;

    }

    if(this->kpoints.size() == 0) 
        throw types::parseError() << types::errinfo_parse("No energy levels found in file.");

    return true;
}


//bool linkWithCubeFiles(std::list<formats::WfnCube> &cubeList){
//    for(Uint k = 0; k < kpoints.size(); ++k){
//        atomistic::EnergyLevels levels = kpoints[k];
//
//        for(Uint level = 1; level <= levels.count(); ++level){
//            Real energy = levels.getLevel(level);
//
//            bool found = false;
//            std::list<formats::WfnCube>::iterator 
//                cubeIt = cubeList.begin(), cubeEnd = cubeList.end();
//            
//            while(cubeIt != cubeEnd){
//
//    Regex bandRegex("b\\d\\d?\\d?");
//
//    Match fermiMatch;
//    regex_search(cubeIt->fileName, fermiMatch, levelRegex);
//    Real fermi = lexical_cast<Real>(fermiMatch);
//
//
//                // Check whether we have the right 
//                if(cubeIt->wfn == level && cubeIt->spin == spin +1){
//                    cubeIt->setEnergy(energy);
//                    newCubeList.push_back(*cubeIt);
//
//                    // exit cube search
//                    found = true;
//                    std::cout << "Found cube file for energy level "
//                                                                                                                                                                                        << level << " at "
//                                                                                                                                                                                                                    << energy << " eV\n";
//                                                                                                                                                                                break;
//



          


}
}
