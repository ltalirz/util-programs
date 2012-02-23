#include "formats.hpp"
#include "io.hpp"
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
bool abs(types::String cubeFile);

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        abs(args["cubefile"].as< types::String >());
    }

    return 0;
}

bool abs(types::String cubeFile) {
    using namespace types;

    formats::Cube c = formats::Cube();
    std::cout << "Reading cube file " << cubeFile << "\n";
    t=clock();
    c.readCubeFile(cubeFile);
    std::cout << "Read cube file in " << (clock() -t)/1000.0 << " ms\n";
    c.grid.abs();
    String outCubeFile = "abs.";
    outCubeFile += io::getFileName(cubeFile);
    std::cout << "Writing absolute value of cube file to " << outCubeFile << "\n";
    c.writeCubeFile(outCubeFile);

    return true;
}


bool parse(int ac, char* av[], po::variables_map& vm) {

    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("version,v", "print version information")
    ("cubefile", po::value<types::String>(), "Gaussian cube file to be scaled.")
    ;

    // Register positional options
    po::positional_options_description p;
    p.add("cubefile", 1);

    // Parse
    po::store(po::command_line_parser(ac,av).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    // Display help message
    if (    !vm.count("cubefile")) {
        std::cout << "Usage: cubeabs [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Feb 22nd 2012\n";
    } else {
        return true;
    }

    return false;
}
