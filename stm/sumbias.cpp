#include "atomistic.h"
#include "types.h"
#include "io.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <string>

#include <boost/program_options.hpp>
#include <boost/spirit/include/qi.hpp>
#include <blitz/array.h>
#include <fftw3.h>


#include <ctime>
time_t t = clock();
//#include <gsl/gsl_const_mksa.h>

namespace po = boost::program_options;
namespace at = atomistic;

// Returns true, if parsing went ok
bool parse(int ac, char* av[], po::variables_map& vm);

// Prepares extrapolation
bool prepare(types::String levelFileName,
             types::String cubeListFileName,
             types::String biasListFileName);

// Returns list with required cubes, if present in cubeList
std::vector<Cube> getRequiredCubes(std::vector<types::Real> requiredLevels,
         std::vector<Cube> cubeList);
// Performs summation
// The cubeList is not passed by reference as the contained cubes
// only have header information. This way, not too many cubes
// have to be kept in memory simultaneously.
bool sum(at::Cube & sum,
         std::vector<types::Real> requiredLevels,
         std::vector<Cube> cubeList);

// Sorts by absolute value
bool absSort (types::Real i, types::Real j) { return ( std::abs(i) < std::abs(j) ); }

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        prepare(args["levels"].as< types::String >(),
                args["cubelist"].as< types::String >(),
                args["biaslist"].as< types::String >()
               );

    }

    return 0;
}

bool prepare(types::String levelFileName,
             types::String cubeListFileName,
             types::String biasListFileName) {
    using namespace types;

    // Read energy levels
    at::Spectrum spectrum = at::Spectrum();
    spectrum.readFromCp2k(levelFileName.c_str());
    spectrum.setFermiZero();
    at::EnergyLevels levels = spectrum.sumSpins();

    // Read descriptions of cube files
    std::ifstream cubeListFile;
    cubeListFile.open(cubeListFileName.c_str());
    if (!cubeListFile.is_open() )
            throw types::fileAccessError() << boost::errinfo_file_name(cubeListFileName);
    
    std::vector<at::Cube> cubeList;
    types::String fileName;
    while (getline(cubeListFile, fileName)) {
        at::Cube t = at::Cube();
        t.readDescription(fileName);
        cubeList.push_back(t);
        
    }
    cubeListFile.close();


    // Read list of bias voltages
    using boost::spirit::double_;
    using boost::spirit::qi::eol;
    using boost::spirit::qi::phrase_parse;
    using boost::spirit::ascii::space;

    types::Binary biasBin;
    io::readBinary(biasListFileName, biasBin);
    typedef types::Binary::const_iterator binIt;
    binIt it = biasBin.begin(), end = biasBin.end();

    std::vector<types::Real> biasList;
    if (! phrase_parse(
        it,
        end,
        double_ % eol,
        space,
        biasList
        )) throw types::parseError() << types::errinfo_parse("list of bias voltages");
     std::sort(biasList.begin(), biasList.end());
     // And separate list into positive and nevative voltages
     std::vector<Real>::const_iterator bIt = biasList.begin(), bEnd = biasList.end();
     std::vector<Real> negList, posList; 
     while(bIt != bEnd){
         if(*bIt > 0) posList.push_back(*bIt);
         else negList.push_back(*bIt);
     }
     

     // Perform summation
     std::vector< std::vector<Real> > biasDomains;
     biasDomains.push_back(posList); biasDomains.push_back(negList);
     std::vector< std::vector<Real> >::iterator listIt = biasDomains.begin(), listEnd = biasDomains.end();
     std::vector<Real> requiredLevels;
     std::vector<Real> levelsCopy = levels;
     while(listIt != listEnd){
         std::vector<Real> biasDomain = *listIt; 
         std::vector<Real>::const_iterator biasIt = biasDomain.begin(), biasEnd = biasDomain.end();
         std::sort(biasIt, biasEnd, absSort);

         at::Cube sum = Cube();
         while(biasIt != biasEnd){
             
             std::vector<Real>::const_iterator levelIt = levelsCopy.begin(), levelEnd = levelsCopy.end();
             while(levelIt != levelEnd){
                 if(*levelIt * *biasIt > 0 && *levelIt * *biasIt <= *biasIt * *biasIt){
                     requiredLevels.push_back(*levelIt);
                     // We don't need to add it again
                     levelsCopy.erase(levelIt);
                 }
                 ++levelIt;
             }
            
             std::vector<Cube> toAdd = get(requiredLevels, cubeList);
             // perform sum of cube files as given by requiredLevels
             sum(sum, toAdd);
             ++biasIt;
         }
         requiredLevels.clear();
         ++listIt;
     }

     return true;
}



bool parse(int ac, char* av[], po::variables_map& vm) {

    // Declare regular options
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "produce help message")
    ("version,v", "print version information")
    ("levels", po::value<types::String>(), "CP2K output file containing energy levels")
    ("cubelist", po::value<types::String>(), "file with list of extrapolated wave function cubes from CP2K")
    ("biaslist", po::value<types::String>(), "file with list of biases that shall be computed")
    ;

    // Register positional options
    po::positional_options_description p;
    p	.add("levels", 1)
    .add("cubelist", 1)
    .add("biaslist", 1)
    ;

    // Parse
    po::store(po::command_line_parser(ac,av).
              options(desc).positional(p).run(), vm);
    po::notify(vm);

    // Display help message
    if (	vm.count("help") ||
            !vm.count("levels") ||
            !vm.count("cubelist") ||
            !vm.count("biaslist")) {
        std::cout << "Usage: sumbias [options]\n";
        std::cout << desc << "\n";
    } else if (vm.count("version")) {
        std::cout << "Feb 1st 2012\n";
    } else {
        return true;
    }

    return false;
}

std::vector<Cube> getRequiredCubes(
        std::vector<types::Real> requiredLevels,
        at::Spectrum spectrum,
        std::vector<at::Cube> cubeList){
    std::vector<Real>::const_iterator levelIt = requiredLevels.begin(), levelEnd = requiredLevels.end();
    std::vector<at::Cube>::const_iterator cubeIt = cubeList.begin(), cubeEnd = cubeList.end();
    std::vector<at::Cube> requiredCubes;
    bool found;
    
    while(levelIt != levelEnd){
        found = false;
        while(cubeIt != cubeEnd){
            if( spectrum[ cubeIt->spin -1] [ cubeIt->wfn -1] == *levelIt ){
                requiredCubes.push_back(*cubeIt);
                cubeList.erase(cubeIt);
                found = true;
                break;
            }
            ++cubeIt;
        }
        if(!found) std::cout << "Required cube file for energy level " << *levelIt << " Ha was not in the list.\n";
        ++levelIt;
    }
 
    return requiredCubes;
}


bool sum(at::Cube & sum,
         std::vector<types::Real> requiredLevels,
         std::vector<Cube> cubeList) {
        
    
//
//    using namespace blitz;
//   
//    std::cout << "--------\n Extrapolating " << cubeFile << "\n";
//
//    t = clock();
//    at::WfnCube wfn = at::WfnCube();
//    wfn.readCubeFile(cubeFile.c_str());
//    std::cout << "Time to read cube : " << (clock() -t)/1000.0 << " ms\n";
//    
//    wfn.energy = spectrum.spins[wfn.spin -1].getLevel(wfn.wfn);
//    types::Uint nX = wfn.grid.directions[0].incrementCount,
//                nY = wfn.grid.directions[1].incrementCount,
//                nZ = endIndex;
//
//    wfn.grid.resize(nX, nY, nZ);
//    Array<types::Real,3> dataArray(
//        &wfn.grid.data[0],
//        shape(nX, nY, nZ),
//        neverDeleteData);
//
//    // Produce z profile
//    std::string zProfile = cubeFile;
//    at::WfnCube wfnSq = wfn;
//    wfnSq.grid.squareValues();
//    zProfile += ".zprofile";
//    wfnSq.writeZProfile(zProfile, "Z profile of sqared wave function\n");
//    std::cout << "Wrote Z profile to " << zProfile << std::endl;
//
//
//    // Get z-plane for interpolation
//    Array<types::Real,2> planeDirect(nX, nY);
//    planeDirect = dataArray(Range::all(), Range::all(), startIndex);
//
//    // Do a real 2 complex fft
//    // Since a(-k)=a(k)^* (or in terms of indices: a[n-k]=a[k]^*)
//    // only a[k], k=0...n/2+1 (division rounded down) are retained
//    types::Uint nXF = nX, nYF = nY/2 + 1;
//    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nXF * nYF);
//    fftw_plan plan_forward = fftw_plan_dft_r2c_2d(nX, nY, planeDirect.data(), out, FFTW_ESTIMATE);
//    fftw_execute(plan_forward);
//
//    Array<types::Complex,2> planeFourier(
//        reinterpret_cast<types::Complex *>(out),
//        shape(nXF, nYF),
//        neverDeleteData);
//
//    // Calculate the exponential prefactors.
//    // Multiplication of Fourier coefficients with these prefactors
//    // propagates them to next z plane.
//    // The prefectors would only need dimensions (nX/2 +1, nY/2 +1)
//    // since they are the same for G and -G. However for the sake of
//    // simplicity of calculation, we prepare them here for (nX, nY/2+1) = (nXF,nYF).
//    t = clock();
//    Array<types::Real,2> prefactors(nXF, nYF);
//    // SI units energy: 2 m E/hbar^2 a0 = 2 E/Ha
//    types::Real energyTerm = 2 * wfn.energy;
//    types::Real deltaZ = wfn.grid.directions[2].incrementVector[2];
//    la::Cell cell = la::Cell(wfn.grid.directions);
//    vector<types::Real> X = cell.vector(0);
//    vector<types::Real> Y = cell.vector(1);
//    types::Real dKX = 2 * M_PI / X[0];
//    types::Real dKY = 2 * M_PI / Y[1];
//
//    // Notice that the order of storage is 0...G -G...-1 for uneven nX
//    // and 0...G-1 G -(G-1)...-1 for even NX
//    prefactors( Range(0, nXF/2), Range::all())= exp(- sqrt(tensor::i * dKX * tensor::i * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);
//    // tensor::i always starts from 0, i.e. it ranges from 0...nXF/2-1 (nX even) or 0...nXF/2 (nX uneven)
//    if(nX % 2 == 1)
//            prefactors( Range(nXF/2 + 1, nXF - 1), Range::all())= exp(- sqrt( (nXF/2 - tensor::i) * dKX * (nXF/2 - tensor::i) * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);
//    else
//            prefactors( Range(nXF/2 + 1, nXF - 1), Range::all())= exp(- sqrt( (nXF/2 -1 - tensor::i) * dKX * (nXF/2 - 1 - tensor::i) * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);
//
//
////    // Produce Fourier coefficient for gnuplot    
////    types::Complex *pt =  reinterpret_cast<types::Complex *>(out);
////    for(types::Uint i = 0; i < nXF; ++i){
////            for(types::Uint j = 0;j< nYF; ++j){
////                    std::cout << std::abs(*pt) << " ";
////                    ++pt;
////            }
////            std::cout << std::endl;
////    }
//
//    // Sequentially update the cube file
//    Array<types::Real,2> tempDirect(nX, nY);
//    Array<types::Complex,2> tempFourier(nXF, nYF);
//    fftw_plan plan_backward;
//    for(int zIndex = startIndex + 1; zIndex <= endIndex; ++zIndex) {
//        // Propagate fourier coefficients
//        planeFourier = prefactors(tensor::i, tensor::j) * planeFourier(tensor::i, tensor::j);
//        
//        // The c2r transform destroys its input
//        tempFourier = planeFourier;
//        plan_backward = fftw_plan_dft_c2r_2d(nX, nY, (fftw_complex*) tempFourier.data(), tempDirect.data(), FFTW_ESTIMATE);
//
//        // Do Fourier-back transform
//        fftw_execute(plan_backward); /* repeat as needed */
//        tempDirect /= nX * nY;
//        
//        // Copy data
//        dataArray(Range::all(), Range::all(), zIndex) = tempDirect(Range::all(), Range::all());
//    }
//
//    fftw_destroy_plan(plan_forward);
//    fftw_destroy_plan(plan_backward);
//    fftw_free(out);
//
//    std::cout << "Time to extrapolate : " << (clock() -t)/1000.0 << " ms\n";
//
//    t = clock();
//    types::String outCubeFile = "extrapolated.";
//    outCubeFile += cubeFile;
//    wfn.writeCubeFile(outCubeFile);
//    std::cout << "Wrote extrapolated cube file to " << outCubeFile << std::endl;
//    std::cout << "Time to write cube : " << (clock() -t)/1000.0 << " ms\n";
//
//    // Produce z profile
//    wfn.grid.squareValues();
//    std::string outZProfile = outCubeFile;
//    outZProfile += ".zprofile";
//    wfn.writeZProfile(outZProfile, "Z profile of squared extrpolated wave function\n");
//    std::cout << "Wrote Z profile of extrapolated cube file to " << outZProfile << std::endl;
//
//    return true;
//}
