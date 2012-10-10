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
             std::vector<types::String> wfnCubes,
             types::String hartreeFileName,
             types::String mode,
             double start,
             double width,
             double isoValue,
             double approachFrom,
             double kPlaneMax,
             types::Uint nLayers,
             double ballRadius
             );

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        prepare(args["levels"].as< types::String >(),
                args["wfncubes"].as< std::vector<types::String> >(),
                args["hartree"].as< types::String >(),
                args["mode"].as< types::String >(),
                args["start"].as< double >(),
                args["width"].as< double >(),
                args["isovalue"].as< double >(),
                args["approach-from"].as< double >(),
                args["k-cutoff"].as< double >(),
                args["nlayers"].as< types::Uint >(),
                args["ball-radius"].as< double >()
               );

    }

    return 0;
}

bool prepare(types::String levelFileName,
             std::vector<types::String> wfnCubes,
             types::String hartreeFileName,
             types::String mode,
             double start,
             double width,
             double isoValue,
             double approachFrom,
             double kPlaneMax,
             types::Uint nLayers,
             double ballRadius
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
    
    std::cout << "Time to prepare : " << (clock() -t)/1000.0 << " ms\n";
    t = clock();
    
    // Iterate over cubes
    stm::WfnExtrapolation::Mode m ;
    types::Real var1;
    if( mode == "plane" ){
        m = stm::WfnExtrapolation::plane;
        var1 = start;
    }
    else if ( mode == "isosurface" ){
        m = stm::WfnExtrapolation::isoSurface;
        var1 = isoValue;
    }
    else if ( mode == "rolling-ball" ){
        m = stm::WfnExtrapolation::rollingBall;
        var1 = ballRadius;
    }
    stm::WfnExtrapolation extrapolation = stm::WfnExtrapolation(
            wfnCubes,
            spectrum,
            hartree,
            m,
            var1,
            width,
            approachFrom,
            kPlaneMax,
            nLayers
            );
    extrapolation.execute();

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
    ("wfncubes", po::value< std::vector<types::String> >(), "List of wavefunction cubes you wish to extrapolate")
    ("hartree", po::value<types::String>(), "cube file of hartree potential from CP2K")
    ("mode", po::value<types::String>()->default_value("plane"), "may be 'plane', 'isosurface', 'rolling-ball'")
    ("start", po::value<double>()->default_value(5), "mode 'plane': distance between extrapolation plane and outermost atom in a.u.")
    ("width", po::value<double>()->default_value(15), "length of extrapolation in a.u.")
    ("isovalue", po::value<double>()->default_value(0.0), "mode 'isosurface': isovalue of the Hartree potential [a.u.] ")
    ("approach-from", po::value<double>()->default_value(-1.0), "mode 'isosurface': z [a.u.] from where you want to go down to find the isosurface (default: top z of cube file).")
//  ("decay-cutoff", po::value<double>()->default_value(2.0), "Maximum decay constant k [1/a.u.] to be retained for z decay  10^(-k*z).")
    ("k-cutoff", po::value<double>()->default_value(3.0), "Restrict basis functions to maximum wave number k = 1/lambda [1/a.u.] in xy-plane.")
    ("nlayers", po::value<types::Uint>()->default_value(1), "Number of layers to fit the wave function values.")
    ("ball-radius", po::value<double>()->default_value(5.0), "mode 'rolling-ball': Radius of ball rolling on top of the atoms [a.u.] ")
    ;

    // Register positional options
    po::positional_options_description p;
    p	.add("input-file", 1) .add("wfncubes", -1);

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
            !vm.count("wfncubes") ||
            !vm.count("width")	) {
        std::cout << "Usage: extrapolate [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "October 9th 2012\n";
    } else if ( vm.count("mode")  && 
                vm["mode"].as< types::String >() != "plane" && 
                vm["mode"].as< types::String >() != "rolling-ball" && 
                vm["mode"].as< types::String >() != "isosurface") {
                std::cout << "Error: invalid mode specified.\n";
    } else if ( (!vm.count("mode")  || vm["mode"].as< types::String >() == "plane") &&
               !(vm.count("start") && vm.count("width"))) {
                std::cout << "Error: Need to specify 'start' and 'width' for mode='plane'.\n";
    } else if ( vm["mode"].as< types::String >() == "isosurface" &&
               !(vm.count("isovalue") && vm.count("width"))) {
                std::cout << "Error: Need to specify 'isovalue' and 'width' for mode='isosurface'.\n";
    } else if ( vm["mode"].as< types::String >() == "rolling-ball" &&
               !(vm.count("ball-radius") && vm.count("width"))) {
                std::cout << "Error: Need to specify 'ball-radius' and 'width' for mode='rolling-ball'.\n";
    } else if ( vm["k-cutoff"].as< types::Real >() < 0 ){
                std::cout << "Error: K-cutoff must be non-negative.\n";
    } else {
        return true;
    }

    return false;
}

