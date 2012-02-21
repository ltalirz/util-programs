#include "formats.h"
#include "io.h"
#include "types.hpp"

#include <vector>
#include <string>

#include <boost/program_options.hpp>

#include <ctime>
time_t t = clock();

namespace po = boost::program_options;

// Returns true, if command line arguments parse ok
bool parse(int ac, char* av[], po::variables_map& vm);

// Returns true if cube file was strided successfully.
bool makeStride(types::String cubeFile,
        std::vector<types::Uint> strides);

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        makeStride(args["cubefile"].as< types::String >(),
                args["stride"].as< std::vector<types::Uint> >());
    }

    return 0;
}

bool makeStride(types::String cubeFile,
        std::vector<types::Uint> strides) {
    using namespace types;

    formats::Cube c = formats::Cube();
    std::cout << "Reading cube file " << cubeFile << "\n";
    t=clock();
    c.readCubeFile(cubeFile);
    std::cout << "Read cube file in " << (clock() -t)/1000.0 << " ms\n";
    c.grid.stride(strides);
    String outCubeFile = "strided.";
    outCubeFile += io::getFileName(cubeFile);
    std::cout << "Writing strided cube file to " << outCubeFile << "\n";
    c.writeCubeFile(outCubeFile);

    return true;
}


bool parse(int ac, char* av[], po::variables_map& vm) {

    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("version,v", "print version information")
    ("cubefile", po::value<types::String>(), "Gaussian cube file to be strided.")
    ("stride", po::value< std::vector<types::Uint> >()->multitoken(), "A sequence of unsigned integers, e.g. 2 2 2.")
    ;

    // Register positional options
    po::positional_options_description p;
    p.add("cubefile", 1).add("stride", -1);

    // Parse
    po::store(po::command_line_parser(ac,av).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    // Display help message
    if (    !vm.count("cubefile") ||
            !vm.count("stride")) {
        std::cout << "Usage: cubestride [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Feb 20th 2012\n";
    } else {
        return true;
    }

    return false;
}
