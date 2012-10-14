#include "atomistic.hpp"
#include "formats.hpp"
#include "types.hpp"
#include "io.hpp"

#include <iostream>
#include <fstream>
#include <iterator>
#include <list>
#include <vector>
#include <algorithm>
#include <string>

#include <boost/format.hpp>
#include <boost/program_options.hpp>

#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi_eol.hpp>

#include <blitz/array.h>
#include <fftw3.h>


#include <ctime>
time_t t = clock();

namespace po = boost::program_options;
namespace at = atomistic;
namespace cp2k = formats::cp2k;
using namespace types;

// Returns true, if command line arguments parse ok
bool parse(int ac, char* av[], po::variables_map& vm);

// Reads and parses the lists suppiled by command line
// and calls sum
bool readLists(types::String levelFileName,
             types::String cubeListFileName,
             types::String biasListFileName,
             bool psisquared);

// Performs summation
// The cubes in cubeList only have header information.
// No more than two cubes are kept in memory simultaneously at any time.
bool sum(std::vector<Real> biasDomain,
        std::list<formats::WfnCube> &cubeList,
        cp2k::Spectrum spectrum,
        bool psiSquared);

// Adds appropriate description and writes cube file
bool write(formats::Cube & sum,
           Real bias);

// Sorts by absolute value
bool absSort(Real a, Real b) { return std::abs(a) < std::abs(b); }




int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        readLists(args["levels"].as< types::String >(),
                args["cubelist"].as< types::String >(),
                args["biaslist"].as< types::String >(),
                args["psisquared"].as< bool >()
               );

    }

    return 0;
}

bool readLists(types::String levelFileName,
        types::String cubeListFileName,
        types::String biasListFileName,
        bool psiSquared) {
    using namespace types;

    // Read energy levels
    cp2k::Spectrum spectrum = cp2k::Spectrum();
    spectrum.readFromCp2k(levelFileName.c_str());
    spectrum.setFermiZero();
    // Want spectrum in eV
    spectrum *= at::units::Ha / at::units::eV;
    std::cout << "Read energy levels from " << levelFileName << "\n";

    // Read descriptions of cube files
    std::ifstream cubeListFile;
    cubeListFile.open(cubeListFileName.c_str());
    if (!cubeListFile.is_open() )
        throw types::fileAccessError() << boost::errinfo_file_name(cubeListFileName);

    std::list<formats::WfnCube> cubeList;
    types::String fileName;
    while (getline(cubeListFile, fileName)) {
        formats::WfnCube t = formats::WfnCube();
        t.readDescription(fileName);
        cubeList.push_back(t);

    }
    cubeListFile.close();
    std::cout << "Read list of cube files from " << cubeListFileName << "\n";


    // Read list of bias voltages
    using boost::spirit::double_;
    using boost::spirit::qi::eol;
    using boost::spirit::qi::phrase_parse;
    using boost::spirit::ascii::space;

    types::Binary biasBin;
    io::readBinary(biasListFileName, biasBin);
    typedef types::Binary::const_iterator binIt;
    binIt it = biasBin.begin(), end = biasBin.end();
    std::vector<Real> biasList;
    if (! phrase_parse(
                it,
                end,
                *double_,
                space,
                biasList
                )) throw types::parseError() << types::errinfo_parse("list of bias voltages");
    std::sort(biasList.begin(), biasList.end());

    // Separate list into positive and nevative voltages
    std::vector<Real>::const_iterator bIt = biasList.begin(), bEnd = biasList.end();
    std::vector<Real> negList, posList; 
    while(bIt != bEnd){
        if(*bIt > 0) posList.push_back(*bIt);
        else negList.push_back(*bIt);
        ++bIt;
    }
    std::vector< std::vector<Real> > biasDomains;
    if(! posList.empty() ) { biasDomains.push_back(posList);}
    if(! negList.empty() ) { biasDomains.push_back(negList);}
    std::cout << "Read list of bias voltages from " << biasListFileName << "\n";


    // Process positive and negative list separately
    std::vector< std::vector<Real> >::iterator listIt = biasDomains.begin(), listEnd = biasDomains.end();
    while(listIt != listEnd){
        std::cout << "Starting summation\n";
        sum(*listIt, cubeList, spectrum, psiSquared);
        ++listIt;

    }

    return true;
}


bool sum(std::vector<Real> biasDomain,
        std::list<formats::WfnCube> &cubeList,
        cp2k::Spectrum spectrum,
        bool psiSquared){

    sort(biasDomain.begin(), biasDomain.end(), absSort);
    std::vector<Real>::iterator biasIt = biasDomain.begin(), 
    biasEnd = biasDomain.end();
    formats::WfnCube sum = formats::WfnCube();

    while(biasIt != biasEnd){
       std::cout << "\nBias " << *biasIt << " V\n--------------\n\n";

        Uint nToSum = 0;
        for(Uint spin = 0; spin < spectrum.spins.size(); ++spin){
            for(Uint level = 1; level <= spectrum.spins[spin].count(); ++level){
                Real energy = spectrum.spins[spin].getLevel(level);

                // If we need this level ...
                if(energy != 1e6 && energy * *biasIt >= 0 && energy * *biasIt <= *biasIt * *biasIt){
                    // Find cube file
                    bool found = false;
                    std::list<formats::WfnCube>::const_iterator
                        cubeIt = cubeList.begin(),
                        cubeEnd = cubeList.end();
                    while(cubeIt != cubeEnd){
                        if(cubeIt->wfn == level && cubeIt->spin == spin +1){
                            // If sum empty, take cube
                            if(sum.grid.data.size() == 0){
                                sum = *cubeIt;
                                sum.readCubeFile();
                                if (!psiSquared) sum.squareValues();
                            }
                            // Else perform summation
                            else{
                                // Make local copy of cube file and then read
                                formats::WfnCube temp = *cubeIt;
                                temp.readCubeFile();
                                if (!psiSquared) temp.squareValues();
                                sum += temp;
                            }
                            
                            // Mark level as used and exit cube search
                            found = true; ++nToSum;
                            spectrum.spins[spin].setLevel(level, 1e6);
                            std::cout << "Added cube file for energy level "
                                << level << " at "
                                << energy << " eV\n";
                            break;
                        }
                        else ++cubeIt;

                    }
                    if(!found) std::cout << "Missing cube file for energy level "
                        << level << " at "
                        << energy << " eV\n";
                }
            }
        }
        if(nToSum == 0) std::cout << "No new cubes for bias " << *biasIt << " V\n";

        write(sum, *biasIt);

        ++biasIt;
    }

    return true;
}


bool parse(int ac, char* av[], po::variables_map& vm) {

    std::string input_file;
    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "produce help message")
    ("version,v", "print version information")
    ("input-file,i", po::value<types::String>(&input_file),"Input file specifying all or several of the following options")
    ("levels", po::value<types::String>(), "CP2K output file containing energy levels")
    ("cubelist", po::value<types::String>(), "file with list of extrapolated wave function cubes from CP2K")
    ("biaslist", po::value<types::String>(), "file with list of biases that shall be computed")
    ("psisquared", po::value<bool>()->default_value(false), "Whether the cube files contain the square of the wave function (and not the wave function itself)")
    ;

    // Register positional options
    po::positional_options_description p;
    p	.add("input-file", -1)
    ;

    // Parse
    po::store(po::command_line_parser(ac,av).
              options(desc).positional(p).run(), vm);
    po::notify(vm);
    
    // If specified, try to parse conig file
    if (vm.count("input-file")){
        std::ifstream ifs(input_file.c_str());
       if(!ifs) throw types::fileAccessError() << boost::errinfo_file_name(input_file);
       store(parse_config_file(ifs, desc), vm);
       notify(vm);
    }

    // Display help message
    if (	vm.count("help") ||
            !vm.count("levels") ||
            !vm.count("cubelist") ||
            !vm.count("biaslist")) {
        std::cout << "Usage: sumbias [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Apr 20th 2012\n";
    } else {
        return true;
    }

    return false;
}

bool write(formats::Cube & sum,
           Real bias) {
    types::String biasString = str(boost::format("%4.3f") % bias);

    types::String description = "Sum of (squared) cube files for STM at bias ";
    description += biasString; 
    description += " V";
    sum.description = description;

    types::String fileName = "bias_";
    fileName += biasString;
    fileName += ".cube";

    std::cout << "Writing file " << fileName << "\n";
    sum.writeCubeFile(fileName);

    return true;
}
