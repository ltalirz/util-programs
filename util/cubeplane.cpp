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
bool calculate(types::String cubeFile, 
               types::Uint dir, 
               types::Uint index);

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        types::Uint dir;
        types::String dirLetter = args["dir"].as< types::String >();
        if( dirLetter == "x") dir = 0;
        else if ( dirLetter == "y") dir = 1;
        else if ( dirLetter == "z") dir = 2;
        
        calculate(args["cubefile"].as< types::String >(),  
                  dir, 
                  args["index"].as< types::Uint >());
    }

    return 0;
}

bool calculate(types::String cubeFile, types::Uint dir, types::Uint index) {
    using namespace types;

    formats::Cube c = formats::Cube();
    std::cout << "Reading cube file " << cubeFile << "\n";
    t=clock();
    c.readCubeFile(cubeFile);
    std::cout << "Read cube file in " << (clock() -t)/1000.0 << " ms\n";
    String outFile = io::getFileName(cubeFile);
    
    std::vector<types::Real> plane;
    c.grid.plane(dir, index, plane);
    outFile += ".plane";
    std::cout << "Writing plane with index " << index << " to " << outFile << "\n";
    std::cout << "Plane dimensions " << plane.size() << std::endl;
    c.grid.writeDirPlane(outFile, plane, dir);

    return true;
}


bool parse(int ac, char* av[], po::variables_map& vm) {

    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("version,v", "print version information")
    ("cubefile", po::value<types::String>(), "Gaussian cube file")
    ("dir", po::value<types::String>(), "=x, y or z")
    ("index", po::value<types::Uint>(), "Index of plane in cube file.")
    ;

    // Register positional options
    po::positional_options_description p;
    p.add("cubefile", 1).add("dir", 2).add("index", 3);

    // Parse
    po::store(po::command_line_parser(ac,av).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    // Display help message
    if (    !vm.count("cubefile")) {
        std::cout << "Usage: cubeplane [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Aug 3rd 2012\n";
    } else if (vm["dir"].as< types::String >() != "x" &&
            vm["dir"].as< types::String >() != "y" &&
            vm["dir"].as< types::String >() != "z") {
        std::cout << "dir must be 'x', 'y' or 'z'";
    } else {
        return true;
    }

    return false;
}
