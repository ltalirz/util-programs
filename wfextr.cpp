#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <boost/program_options.hpp>

#include "atomistic.h"

namespace po = boost::program_options;
namespace at = atomistic;

// Returns true, if parsing went ok
bool parse(int ac, char* av[], po::variables_map& vm);

// Handles extrapolation
bool extrapolate(std::string levelFile,
                 std::string cubeListFile,
                 std::string hartreeFile,
                 double start,
                 double width);

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        extrapolate(args["levels"].as< std::string >(),
                    args["cubelist"].as< std::string >(),
                    args["hartree"].as< std::string >(),
                    args["start"].as< double >(),
                    args["width"].as< double >()
                   );

    }

    return 0;
}

bool extrapolate(std::string levelFile,
                 std::string cubeListFile,
                 std::string hartreeFile,
                 double start,
                 double width) {

    // Read energy levels
    at::Spectrum spectrum = at::Spectrum();
    spectrum.readFromCp2k(levelFile.c_str());

    // Read hartree cube file
    at::Cube cube = at::Cube();
    cube.readCubeFile(hartreeFile.c_str());
    cube.grid.squareValues();
    std::vector<double> line;
    cube.grid.sumXY(line);
    std::vector<double>::iterator it;
//    for(it = line.begin(); it!=line.end(); ++it) {
//        std::cout << *it << std::endl;
//    }


    return true;
}


bool parse(int ac, char* av[], po::variables_map& vm) {

    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "produce help message")
    ("version,v", "print version information")
    ("levels", po::value<std::string>(), "CP2K output file containing energy levels")
    ("cubelist", po::value<std::string>(), "file with list of wavefunction cubes you wish to extrapolate")
    ("hartree", po::value<std::string>(), "cube file of hartree potential from CP2K")
    ("start", po::value<double>(), "distance between extrapolation plane and outermost atom in a.u.")
    ("width", po::value<double>(), "length of extrapolation in a.u.")
    ;

    // Register positional options
    po::positional_options_description p;
    p	.add("levels", 1)
    .add("cubelist", 1)
    .add("hartree", 1)
    .add("start", 1)
    .add("width", 1);

    // Parse
    po::store(po::command_line_parser(ac,av).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    // Display help message
    if (	vm.count("help") ||
            !vm.count("levels") ||
            !vm.count("hartree") ||
            !vm.count("start") ||
            !vm.count("width")	) {
        std::cout << "Usage: wfextr.x [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Dez 5th 2011\n";
    } else {
        return true;
    }

    return false;
}

