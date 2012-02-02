#include "atomistic.h"
#include "types.h"
#include "io.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <string>

#include <boost/program_options.hpp>
#include <boost/spirit/include/qi.hpp>
#include <blitz/array.h>
#include <fftw3.h>


#include <ctime>
time_t t = clock();
//#include <gsl/gsl_const_mksa.h>

namespace po = boost::program_options;
namespace at = atomistic;

// Returns true, if parsing went ok
bool parse(int ac, char* av[], po::variables_map& vm);

// Prepares extrapolation
bool prepare(types::String levelFileName,
             types::String cubeListFileName,
             types::String biasListFileName);

// Returns list with required cubes, if present in cubeList
std::vector<at::WfnCube> getRequiredCubes(
        std::vector<types::Real> requiredLevels,
        at::Spectrum spectrum,
        std::vector<at::WfnCube> cubeList);
// Performs summation
// The cubeList is not passed by reference as the contained cubes
// only have header information. This way, not too many cubes
// have to be kept in memory simultaneously.
bool sum(at::Cube & sum,
         std::vector<at::WfnCube> toAdd);

bool write(at::Cube & sum,
           types::Real bias);

// Sorts by absolute value
bool absSort (types::Real i, types::Real j) { return ( std::abs(i) < std::abs(j) ); }

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        prepare(args["levels"].as< types::String >(),
                args["cubelist"].as< types::String >(),
                args["biaslist"].as< types::String >()
               );

    }

    return 0;
}

bool prepare(types::String levelFileName,
             types::String cubeListFileName,
             types::String biasListFileName) {
    using namespace types;

    // Read energy levels
    at::Spectrum spectrum = at::Spectrum();
    spectrum.readFromCp2k(levelFileName.c_str());
    spectrum.setFermiZero();
    at::EnergyLevels levels = spectrum.sumSpins();

    // Read descriptions of cube files
    std::ifstream cubeListFile;
    cubeListFile.open(cubeListFileName.c_str());
    if (!cubeListFile.is_open() )
            throw types::fileAccessError() << boost::errinfo_file_name(cubeListFileName);
    
    std::vector<at::WfnCube> cubeList;
    types::String fileName;
    while (getline(cubeListFile, fileName)) {
        at::WfnCube t = at::WfnCube();
        t.readDescription(fileName);
        cubeList.push_back(t);
        
    }
    cubeListFile.close();


    // Read list of bias voltages
    using boost::spirit::double_;
    using boost::spirit::qi::eol;
    using boost::spirit::qi::phrase_parse;
    using boost::spirit::ascii::space;

    types::Binary biasBin;
    io::readBinary(biasListFileName, biasBin);
    typedef types::Binary::const_iterator binIt;
    binIt it = biasBin.begin(), end = biasBin.end();

    std::vector<types::Real> biasList;
    if (! phrase_parse(
        it,
        end,
        double_ % eol,
        space,
        biasList
        )) throw types::parseError() << types::errinfo_parse("list of bias voltages");
     std::sort(biasList.begin(), biasList.end());
     // And separate list into positive and nevative voltages
     std::vector<Real>::const_iterator bIt = biasList.begin(), bEnd = biasList.end();
     std::vector<Real> negList, posList; 
     while(bIt != bEnd){
         if(*bIt > 0) posList.push_back(*bIt);
         else negList.push_back(*bIt);
     }
     

     // Perform summation
     std::vector< std::vector<Real> > biasDomains;
     biasDomains.push_back(posList); biasDomains.push_back(negList);
     std::vector< std::vector<Real> >::iterator listIt = biasDomains.begin(), listEnd = biasDomains.end();
     std::vector<Real> requiredLevels;
     std::vector<Real> levelsCopy = levels.levels;
     while(listIt != listEnd){
         std::vector<Real> biasDomain = *listIt; 
         std::vector<Real>::iterator biasIt = biasDomain.begin(), biasEnd = biasDomain.end();
         std::sort(biasIt, biasEnd, absSort);

         at::Cube sum = at::Cube();
         while(biasIt != biasEnd){
             
             std::vector<Real>::iterator levelIt = levelsCopy.begin(), levelEnd = levelsCopy.end();
             while(levelIt != levelEnd){
                 if(*levelIt * *biasIt > 0 && *levelIt * *biasIt <= *biasIt * *biasIt){
                     requiredLevels.push_back(*levelIt);
                     // We don't need to add it again
                     levelsCopy.erase(levelIt);
                 }
                 ++levelIt;
             }
            
             std::vector<at::WfnCube> toAdd = getRequiredCubes(requiredLevels, spectrum, cubeList);
             // perform sum of cube files as given by requiredLevels
             ::sum(sum, toAdd);
             write(sum, *biasIt);
             ++biasIt;
         }
         requiredLevels.clear();
         ++listIt;
     }

     return true;
}



bool parse(int ac, char* av[], po::variables_map& vm) {

    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "produce help message")
    ("version,v", "print version information")
    ("levels", po::value<types::String>(), "CP2K output file containing energy levels")
    ("cubelist", po::value<types::String>(), "file with list of extrapolated wave function cubes from CP2K")
    ("biaslist", po::value<types::String>(), "file with list of biases that shall be computed")
    ;

    // Register positional options
    po::positional_options_description p;
    p	.add("levels", 1)
    .add("cubelist", 1)
    .add("biaslist", 1)
    ;

    // Parse
    po::store(po::command_line_parser(ac,av).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    // Display help message
    if (	vm.count("help") ||
            !vm.count("levels") ||
            !vm.count("cubelist") ||
            !vm.count("biaslist")) {
        std::cout << "Usage: sumbias [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Feb 1st 2012\n";
    } else {
        return true;
    }

    return false;
}

// Both required levels and spectrum have fermi at zero
std::vector<at::WfnCube> getRequiredCubes(
        std::vector<types::Real> requiredLevels,
        at::Spectrum spectrum,
        std::vector<at::WfnCube> cubeList){

    using namespace types;
    std::vector<Real>::const_iterator levelIt = requiredLevels.begin(), levelEnd = requiredLevels.end();
    std::vector<at::WfnCube>::iterator cubeIt = cubeList.begin(), cubeEnd = cubeList.end();
    std::vector<at::WfnCube> requiredCubes;
    bool found;

    while(levelIt != levelEnd){
        found = false;
        while(cubeIt != cubeEnd){
            if( spectrum.spins[ cubeIt->spin -1].getLevel(cubeIt->wfn) == *levelIt ){
                requiredCubes.push_back(*cubeIt);
                cubeList.erase(cubeIt);
                found = true;
                break;
            }
            ++cubeIt;
        }
        if(!found) std::cout << "Required cube file for energy level " << *levelIt << " Ha was not in the list.\n";
        ++levelIt;
    }
 
    return requiredCubes;
}


bool sum(at::Cube & sum,
         std::vector<at::WfnCube> toAdd) {
    
    std::vector<at::WfnCube>::const_iterator cubeIt = toAdd.begin(), cubeEnd = toAdd.end();
    while(cubeIt != cubeEnd){
        at::Cube tempCube = *cubeIt;
        tempCube.readCubeFile();
        sum += tempCube;
    }

    return true;
} 

bool write(at::Cube & sum,
           types::Real bias) {
    types::String title = "Sum of cube files for STM at bias ";
    title += bias;
    title += " V\n";
    sum.title.insert(sum.title.begin(), title.begin(), title.end());

    types::String fileName = "bias_";
    fileName += bias;
    fileName += ".cube";
    sum.writeCubeFile(fileName);

    return true;
}
