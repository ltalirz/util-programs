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
        bool psisquared,
        std::vector<Real> values,
        stm::StmCube::Mode m
        );

/***********************************************
  Implementations
***********************************************/

const int CONSTANT_CURRENT = 0;
const int CONSTANT_HEIGHT =  1;


int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        std::vector<Real> v;
        stm::StmCube::Mode m;
    
        // Iterate over cubes
        if(args.count("isovalues")){
                v = args["isovalues"].as< std::vector<Real> >();
                m = stm::StmCube::CONSTANT_CURRENT;
         }
         else{
                v = args["zvalues"].as< std::vector<Real> >();
                m = stm::StmCube::CONSTANT_HEIGHT;
         }
         doStm(
                args["cubelist"].as< types::String >(),
                args["psisquared"].as< bool >(),
                v,
                m
              );
    }

    return 0;
}


bool doStm(types::String cubeListFileName, 
           bool psisquared,
           std::vector<Real> values,
           stm::StmCube::Mode m
           ){
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

        std::vector<Real>::const_iterator isoIt = values.begin(),
            isoEnd = values.end();
        while(isoIt != isoEnd){
            types::String igorFileName;
            if (m == stm::StmCube::CONSTANT_CURRENT){
                temp.setIsoValue(*isoIt);
                igorFileName = io::getFileName(fileName);
                igorFileName += str(boost::format(".%2.1e.igor") % *isoIt);
            }
            else{
                temp.setZValue(*isoIt);
                igorFileName = io::getFileName(fileName);
                igorFileName += str(boost::format(".%1.2d.igor") % *isoIt);
            }
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
    ("isovalues", po::value< std::vector<Real> >()->multitoken(), "The isovalues defining the density isosurface of an STM image in constantr-current mode. 1E-7 is typically a good start.")
    ("zvalues", po::value< std::vector<Real> >()->multitoken(), "The height above the topmost atoms for an STM image in constant-height mode. 8 a.u. is typically reasonable.")
    ("psisquared", po::value<bool>()->default_value(false), "Whether the cube files contain the square of the wave function (and not the wave function itself)")
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
            !vm.count("isovalues") && !vm.count("zvalues") ||
            !vm.count("cubelist")
       ){
        std::cout << "Usage: stm [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "September 26th 2012\n";
    } else {
        return true;
    }

    return false;
}
