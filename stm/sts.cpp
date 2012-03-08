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
bool readLists(types::String levelFileName,
             types::String cubeListFileName,
             types::String out,
             types::Real eMin,
             types::Real eMax,
             types::Real deltaE,
             types::Real FWHM,
             stm::StsCube::ModeFlag modeFlag,
             types::String modeParameter);

// Inserts info about energy levels in cube file descriptions
bool insertEnergyLevels( 
        std::list<formats::WfnCube> &cubeList,
        formats::cp2k::Spectrum spectrum,
        types::Real eMin,
        types::Real eMax,
        types::Real sigma);

/***********************************************
  Implementations
***********************************************/

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {

        // Process mode-related input
        stm::StsCube::ModeFlag modeFlag;
        String modeParameter;
        
        if (args["mode"].as< String >() == "constant_z"){
            modeFlag = stm::StsCube::CONSTANT_Z;
            modeParameter = args["height"].as< String >() ;
        } else if (args["mode"].as< String >() == "profile"){
            modeFlag = stm::StsCube::PROFILE;
            modeParameter = args["zprofile"].as< String >() ;
        }
        
        readLists(args["levels"].as< types::String >(),
                args["cubelist"].as< types::String >(),
                args["out"].as< types::String >(),
                args["emin"].as< types::Real >(),
                args["emax"].as< types::Real >(),
                args["delta-e"].as< types::Real >(),
                args["FWHM"].as< types::Real >(),
                modeFlag,
                modeParameter
                );

    }

    return 0;
}

bool readLists(types::String levelFileName,
             types::String cubeListFileName,
             types::String outFileName,
             types::Real eMin,
             types::Real eMax,
             types::Real deltaE,
             types::Real FWHM,
             stm::StsCube::ModeFlag modeFlag,
             types::String modeParameter){
    using namespace types;

    // Read energy levels
    std::cout << "Reading energy levels from " << levelFileName << "\n";
    formats::cp2k::Spectrum spectrum = formats::cp2k::Spectrum();
    spectrum.readFromCp2k(levelFileName.c_str());
    spectrum.setFermiZero();
    // Want spectrum in eV
    spectrum *= at::units::Ha / at::units::eV;

    // Read descriptions of cube files
    std::cout << "Reading list of cube files from " << cubeListFileName << "\n";
    std::ifstream cubeListFile;
    cubeListFile.open(cubeListFileName.c_str());
    if (!cubeListFile.is_open() )
        throw types::fileAccessError() << boost::errinfo_file_name(cubeListFileName);

    std::list<formats::WfnCube> cubeList;
    types::String fileName;
    while (getline(cubeListFile, fileName)) {
        formats::WfnCube t = formats::WfnCube();
        t.readDescription(fileName);
        cubeList.push_back(t);

    }
    cubeListFile.close();
    insertEnergyLevels(cubeList, spectrum, eMin, eMax, FWHM);

    std::cout << "Done with preparation\n------------\n";
    stm::StsCube mySts = stm::StsCube(
        cubeList,
        eMin,
        eMax,
        deltaE,
        FWHM,
        modeFlag,
        modeParameter);
    std::cout << "Writing STS cube file " << outFileName << "\n";
    mySts.writeCubeFile(outFileName);
    
    return true;
}


bool insertEnergyLevels( 
        std::list<formats::WfnCube> &cubeList,
        formats::cp2k::Spectrum spectrum,
        types::Real eMin,
        types::Real eMax,
        types::Real FWHM){

    std::list<formats::WfnCube> newCubeList;
    // We also take levels that are sigma*3 away from the border
    types::Real delta = 3 * FWHM;

    for(Uint spin = 0; spin < spectrum.spins.size(); ++spin){
        at::EnergyLevels levels = spectrum.spins[spin];
        
        for(Uint level = 1; level <= levels.count(); ++level){
            Real energy = levels.getLevel(level);

            // If we need this level ...
            if( energy >= eMin - delta && energy <= eMax + delta){

                // Find cube file
                bool found = false;
                std::list<formats::WfnCube>::iterator
                    cubeIt = cubeList.begin(),
                           cubeEnd = cubeList.end();
                while(cubeIt != cubeEnd){
                    if(cubeIt->wfn == level && cubeIt->spin == spin +1){
                        cubeIt->setEnergy(energy);
                        newCubeList.push_back(*cubeIt);

                        // exit cube search
                        found = true;
                        std::cout << "Found cube file for energy level "
                            << level << " at "
                            << energy << " eV\n";
                        break;
                    }
                    else ++cubeIt;

                }
                if(!found) std::cout << "Missing cube file for energy level "
                    << level << " at "
                        << energy << " eV\n";
            }
        }
    }
    
    cubeList = newCubeList;

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
    ("levels", po::value<types::String>(), "CP2K output file containing energy levels")
    ("cubelist", po::value<types::String>(), "file with list of extrapolated wave function cubes from CP2K")
    ("out", po::value<types::String>()->default_value("sts.cube"), "filename of sts cube file")
    ("emin", po::value<types::Real>(), "Minimum bias for STS")
    ("emax", po::value<types::Real>(), "Maximum bias for STS")
    ("delta-e", po::value<types::Real>()->default_value(0.01), "Bias step for STS [V]")
    ("FWHM", po::value<types::Real>()->default_value(0.2), "FWHM of Gaussian broadening [V]. FWHM = 2.355 sigma.")
    ("mode", po::value<types::String>()->default_value("constant_z"), "Scanning mode, may be 'constant_z' or 'read_profile'")
    ("height", po::value<types::String>(), "Height [a0] above top surface, where STS shall be performed (mode = 'constant_z')")
    ("zprofile", po::value<types::String>(), ".igor file containing the z profile (a.u.) on which to perform the STS, e.g. STM profile (mode = 'profile')")
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
            !vm.count("levels") ||
            !vm.count("cubelist") ||
            !vm.count("emax") ||
            !vm.count("emin") ||
            !vm.count("delta-e") ||
            !vm.count("FWHM") ||
            !vm.count("mode")
       ){
        std::cout << "Usage: sts [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Mar 1st 2012\n";
    } else if (vm["mode"].as< String >() == "constant_z" && !vm.count("height")) {
                std::cout << "Error: Need to specify 'height' for mode='constant_z'.\n";
    } else if (vm["mode"].as< String >() == "profile" && !vm.count("zprofile")) {
                std::cout << "Error: Need to specify 'zprofile' for mode='profile'.\n";
    } else if (vm["emin"].as< Real >() > vm["emax"].as< Real >()){
        std::cout << "Error: emin > emax\n";
    } else {
        return true;
    }

    return false;
}
