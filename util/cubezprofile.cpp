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
bool zprofile(types::String cubeFile);

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        zprofile(args["cubefile"].as< types::String >());
    }

    return 0;
}

bool zprofile(types::String cubeFile) {
    using namespace types;

    formats::Cube c = formats::Cube();
    std::cout << "Reading cube file " << cubeFile << "\n";
    t=clock();
    c.readCubeFile(cubeFile);
    std::cout << "Read cube file in " << (clock() -t)/1000.0 << " ms\n";
    //c.grid.squareValues();
    String outFile = io::getFileName(cubeFile);
    outFile += ".zprofile";
    std::cout << "Writing z profile to " << outFile << "\n";
    c.writeZProfile(outFile);

    return true;
}


bool parse(int ac, char* av[], po::variables_map& vm) {

    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("version,v", "print version information")
    ("cubefile", po::value<types::String>(), "Gaussian cube file")
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
        std::cout << "Usage: cubezprofile [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Mar 19th 2012\n";
    } else {
        return true;
    }

    return false;
}
