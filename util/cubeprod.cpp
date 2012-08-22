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
bool prod(
        types::String cube1,
        std::vector<types::String> cube2
        );

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        prod(args["cube1"].as< types::String >(),
             args["cube2"].as< std::vector<types::String> >()
             );
    }

    return 0;
}

bool prod(
        types::String cube1,
        std::vector<types::String> cube2
        ){
    using namespace types;

    formats::Cube c1 = formats::Cube(), c2 = formats::Cube();

    std::cout << "Reading cube file " << cube1 << "\n";
    t=clock();
    c1.readCubeFile(cube1);
    std::cout << "Read cube file in " << (clock() -t)/1000.0 << " ms\n";

    std::vector<String>::const_iterator c2it = cube2.begin(), c2end = cube2.end();
    while(c2it != c2end){

        std::cout << "Reading cube file " << *c2it << "\n";
        t=clock();
        c2.readCubeFile(*c2it);
        std::cout << "Read cube file in " << (clock() -t)/1000.0 << " ms\n";

        std::cout << cube1 << " times " << *c2it << " = " 
                  << c1.grid * c2.grid *c1.grid.volumeElement() << "\n";
        ++c2it;
    }

    return true;
}


bool parse(int ac, char* av[], po::variables_map& vm) {

    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("version,v", "print version information")
    ("cube1", po::value<types::String>(), "Gaussian cube file.")
    ("cube2", po::value< std::vector<types::String> >(), "Gaussian cube file(s). ")
    ;

    // Register positional options
    po::positional_options_description p;
    p.add("cube1", 1).add("cube2", -1);

    // Parse
    po::store(po::command_line_parser(ac,av).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    // Display help message
    if (    !vm.count("cube1") ||
            !vm.count("cube2")) {
        std::cout << "Usage: cubeprod [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Aug 22nd 2012\n";
    } else {
        return true;
    }

    return false;
}
