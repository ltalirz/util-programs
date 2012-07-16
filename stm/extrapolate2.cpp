#include "atomistic.hpp"
#include "formats.hpp"
#include "types.hpp"
#include "io.hpp"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>

#include <boost/program_options.hpp>
#include <blitz/array.h>


#include <ctime>
time_t t = clock();
//#include <gsl/gsl_const_mksa.h>

namespace po = boost::program_options;
namespace at = atomistic;
namespace cp2k = formats::cp2k;
namespace stm = formats::stm;

// Returns true, if parsing went ok
bool parse(int ac, char* av[], po::variables_map& vm);

// Prepares extrapolation
bool prepare(types::String levelFileName,
             types::String cubeListFileName,
             types::String hartreeFileName,
             types::String mode,
             double start,
             double width,
             double isoValue
             );

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        prepare(args["levels"].as< types::String >(),
                args["cubelist"].as< types::String >(),
                args["hartree"].as< types::String >(),
                args["mode"].as< types::String >(),
                args["start"].as< double >(),
                args["width"].as< double >(),
                args["isoValue"].as< double >()
               );

    }

    return 0;
}

bool prepare(types::String levelFileName,
             types::String cubeListFileName,
             types::String hartreeFileName,
             types::String mode,
             double start,
             double width,
             double isoValue
            ){

    // Read energy levels
    cp2k::Spectrum spectrum = cp2k::Spectrum();
    spectrum.readFromCp2k(levelFileName.c_str());

    // Read hartree cube file
    formats::Cube hartree = formats::Cube();
    hartree.readCubeFile(hartreeFileName.c_str());

    // Write Zprofile of hartree potential
    std::string hartreeZProfile = io::getFileName(hartreeFileName);
    hartreeZProfile += ".zprofile";
    std::cout << "Writing Z profile of Hartree potential to " 
        << hartreeZProfile << std::endl;
    hartree.writeZProfile(hartreeZProfile);
    
    // Read file with list of cubes
    std::ifstream cubeListFile;
    cubeListFile.open(cubeListFileName.c_str());
    if (!cubeListFile.is_open() )
            throw types::fileAccessError() 
                << boost::errinfo_file_name(cubeListFileName);
    
    std::vector<types::String> cubeList;
    types::String fileName;
    while (getline(cubeListFile, fileName)) {
        cubeList.push_back(fileName);
    }
    cubeListFile.close();

    std::cout << "Time to prepare : " << (clock() -t)/1000.0 << " ms\n";
    t = clock();
    
    // Iterate over cubes
    stm::WfnExtrapolation::Mode m = ( mode == "constantZ" ) 
                                    ? stm::WfnExtrapolation::constantZ 
                                    : stm::WfnExtrapolation::isoValue;
    types::Real var1 = ( m == constantZ ) ? start : isoValue;
    for(std::vector<types::String>::const_iterator it = cubeList.begin();
            it != cubeList.end(); ++it) {
        stm::WfnExtrapolation extrapolation = WfnExtrapolation(
                *it,
                spectrum,
                hartree,
                m,
                var1,
                width);
        extrapolation.execute();
    }

    return true;
}




bool parse(int ac, char* av[], po::variables_map& vm) {
    types::String input_file;
    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "produce help message")
    ("version,v", "print version information")
    ("input-file,i", po::value<types::String>(&input_file), "Input file specifying all or several of the following options")
     ("levels", po::value<types::String>(), "CP2K output file containing energy levels")
    ("cubelist", po::value<types::String>(), "file with list of wavefunction cubes you wish to extrapolate")
    ("hartree", po::value<types::String>(), "cube file of hartree potential from CP2K")
    ("mode", po::value<types::String>(), "may be 'constantZ' or 'isoSurface'")
    ("start", po::value<double>(), "distance between extrapolation plane and outermost atom in a.u.")
    ("width", po::value<double>(), "length of extrapolation in a.u.")
    ("isoValue", po::value<double>(), "isoValue of wavefunction [a.u.] for mode 'isoSurface'")
    ;

    // Register positional options
    po::positional_options_description p;
    p	.add("input-file", -1);

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
            !vm.count("hartree") ||
            !vm.count("width")	) {
        std::cout << "Usage: extrapolate [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "July 15th 2012\n";
    } else if ( vm.count("mode") 
            && (vm["mode"].as< String >() != "constantZ" || vm["mode"].as< String >() != "isoSurface")) {
                std::cout << "Error: invalid mode specified.\n";
    } else if ( (!vm.count("mode")  || vm["mode"].as< String >() == "constantZ") &&
               !(vm.count("start") && vm.count("width"))) {
                std::cout << "Error: Need to specify 'start' and 'width' for mode='constantZ'.\n";
    } else if ( vm["mode"].as< String >() == "isoSurface" &&
               !(vm.count("isoValue") && vm.count("width"))) {
                std::cout << "Error: Need to specify 'isoValue' and 'width' for mode='isoSurface'.\n";
    } else {
        return true;
    }

    return false;
}

