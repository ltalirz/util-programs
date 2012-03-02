#include "atomistic.hpp"
#include "formats.hpp"
#include "types.hpp"
#include "io.hpp"

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <string>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include <ctime>
time_t t = clock();

namespace po = boost::program_options;
namespace at = atomistic;
namespace stm = formats::stm;
using namespace types;

/***********************************************
  Declarations
***********************************************/

// Returns true, if command line arguments parse ok
bool parse(int ac, char* av[], po::variables_map& vm);

// Reads and parses the lists suppiled by command line
// and calls sum
bool doStm(
        types::String cubeListFileName,
        std::vector<Real> isoValues
        );

/***********************************************
  Implementations
***********************************************/

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        doStm(
                args["cubelist"].as< types::String >(),
                args["isovalues"].as< std::vector<Real> >()
           );
    }

    return 0;
}


bool doStm(types::String cubeListFileName, std::vector<Real> isoValues){
    using namespace types;

    // Read cube files
    std::cout << "Reading list of cube files from " << cubeListFileName << "\n";
    std::ifstream cubeListFile;
    cubeListFile.open(cubeListFileName.c_str());
    if (!cubeListFile.is_open() )
        throw types::fileAccessError() << boost::errinfo_file_name(cubeListFileName);

    types::String fileName;

    while(getline(cubeListFile, fileName)) {
        stm::StmCube temp = stm::StmCube();
        std::cout << "Reading cube file " << fileName << "\n";
        temp.readCubeFile(fileName);

        std::vector<Real>::const_iterator isoIt = isoValues.begin(),
            isoEnd = isoValues.end();
        while(isoIt != isoEnd){
            temp.setIsoLevel(*isoIt);
            types::String igorFileName = io::getFileName(fileName);
            igorFileName += str(boost::format(".%2.1e.igor") % *isoIt);
            std::cout << "Writing " << igorFileName << "\n";
            temp.writeIgorFile(igorFileName);
            ++isoIt;
        } 

    }
    cubeListFile.close();

    
    return true;
}


bool parse(int ac, char* av[], po::variables_map& vm) {

    std::string input_file;
    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "produce help message")
    ("version,v", "print version information")
    ("input-file,i", po::value<types::String>(&input_file), "Input file specifying all or several of the following options")
    ("cubelist", po::value<types::String>(), "file with list of extrapolated wave function cubes from CP2K")
    ("isovalues", po::value< std::vector<Real> >()->multitoken(), "The isovalues for the STM image. 1E-7 is typically a good start.")
    ;

    // Register positional options
    po::positional_options_description p;
    p.add("input-file",-1);

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

    // Check if all necessary things are specified
    if (	vm.count("help") ||
            !vm.count("cubelist") ||
            !vm.count("isovalues")
       ){
        std::cout << "Usage: sts [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Mar 1st 2012\n";
    } else {
        return true;
    }

    return false;
}
