#include "atomistic.h"
#include "types.h"
#include "io.h"

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string>

#include <boost/program_options.hpp>
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
             types::String hartreeFileName,
             double start,
             double width);
// Handles extrapolation
/**
 * Do a 2d Fourier transform at a given z-plane,
 * then propagate the Fourier coefficients \f$a_G\f$ according to
 * \f$ a_G(z) = a_G(z_0) e^{-\lambda (z-z_0)}\f$ with
 * \f$\lambda = \sqrt{G^2-\frac{2m}{\hbar^2}(E_n-V_{vac})}\f$
 */
bool extrapolate(types::String cubeFile,
                 types::Uint zStartIndex,
                 types::Uint zEndIndex,
                 const at::Spectrum &spectrum);

int main(int ac, char* av[]) {

    po::variables_map args;
    if (parse(ac, av, args)) {
        prepare(args["levels"].as< types::String >(),
                args["cubelist"].as< types::String >(),
                args["hartree"].as< types::String >(),
                args["start"].as< double >(),
                args["width"].as< double >()
               );

    }

    return 0;
}

bool prepare(types::String levelFileName,
             types::String cubeListFileName,
             types::String hartreeFileName,
             types::Real start,
             types::Real width) {

    // Read energy levels
    at::Spectrum spectrum = at::Spectrum();
    spectrum.readFromCp2k(levelFileName.c_str());

    // Read hartree cube file and shift levels by vacuum level
    at::Cube hartree = at::Cube();
    hartree.readCubeFile(hartreeFileName.c_str());
    types::Real vacuumLevel = hartree.grid.getDataPoint(0,0,hartree.grid.directions[2].incrementCount);
    spectrum.shift(-vacuumLevel);

    std::string hartreeZProfile = hartreeFileName;
    hartreeZProfile += ".zprofile";
    hartree.writeZProfile(hartreeZProfile);

    // Find highest z-coordinate
    // Note: The slowest index in cube file format is x, so in terms of storage
    //       modifying the number of x-coordinates is easiest.
    //       Sadly, physicists like z.
    //       One could reorder the cube data, such that z=>slowest, then do the
    //       extrapolation and reorder again.
    std::vector< at::Atom >::const_iterator it = hartree.atoms.begin(), end = hartree.atoms.end();
    types::Real zTop = it->coordinates[2];
    ++it;
    while(it != end) {
        if( it->coordinates[2] > zTop ) {
            zTop = it->coordinates[2];
        }
//	 else if ( *it.coordinates[2] < zMin ){
//		 zMin = *it.coordinates[2];
//	 }
        ++it;
    }
    // Might do some testing about whether there is enough space to expand
    // for width

    // Define region of interpolation
    std::vector<types::Real> tempvector;
    tempvector.push_back(0);
    tempvector.push_back(0);
    tempvector.push_back(zTop + start);
    std::vector<types::Uint> indices;
    hartree.grid.getNearestIndices(tempvector, indices);
    types::Uint zStartIndex = indices[2];
    tempvector[2] += width;
    hartree.grid.getNearestIndices(tempvector, indices);
    types::Uint zEndIndex = indices[2];

    // Read file with list of cubes
    std::ifstream cubeListFile;
    cubeListFile.open(cubeListFileName.c_str());
    if (!cubeListFile.is_open() )
            throw types::fileAccessError() << boost::errinfo_file_name(cubeListFileName);
    
    std::vector<types::String> cubeList;
    types::String fileName;

    while (getline(cubeListFile, fileName)) {
        cubeList.push_back(fileName);
    }
    cubeListFile.close();

    std::cout << "Time to prepare : " << (clock() -t)/1000.0 << " ms\n";
    t = clock();
    
    // Iterate over cubes
    for(std::vector<types::String>::const_iterator it = cubeList.begin();
            it != cubeList.end(); ++it) {
        extrapolate(*it,zStartIndex,zEndIndex,spectrum);
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
    ("cubelist", po::value<types::String>(), "file with list of wavefunction cubes you wish to extrapolate")
    ("hartree", po::value<types::String>(), "cube file of hartree potential from CP2K")
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
        std::cout << "Jan 19th 2012\n";
    } else {
        return true;
    }

    return false;
}


bool extrapolate(
    types::String cubeFile,
    types::Uint startIndex,
    types::Uint endIndex,
    const at::Spectrum &spectrum) {

    using namespace blitz;
    
    t = clock();
    at::WfnCube wfn = at::WfnCube();
    wfn.readCubeFile(cubeFile.c_str());
    std::cout << "Time to read cube : " << (clock() -t)/1000.0 << " ms\n";
    
    wfn.energy = spectrum.spins[wfn.spin -1].getLevel(wfn.wfn);
    types::Uint nX = wfn.grid.directions[0].incrementCount,
                nY = wfn.grid.directions[1].incrementCount,
                nZ = endIndex;

    
    wfn.grid.resize(nX, nY, nZ);
    Array<types::Real,3> dataArray(
        &wfn.grid.data[0],
        shape(nX, nY, nZ),
        neverDeleteData);

    // Produce z profile
    std::string zProfile = cubeFile;
    at::WfnCube wfnSq = wfn;
    wfnSq.grid.squareValues();
    zProfile += ".zprofile";
    wfnSq.writeZProfile(zProfile, "Z profile of sqared wave function\n");


    // Get z-plane for interpolation
    Array<types::Real,2> planeDirect(dataArray(Range::all(), Range::all(), startIndex));

    // Do a real 2 complex fft
    // Since a(-k)=a(k)^* (or in terms of indices: a[n-k]=a[k]^*)
    // only a[k], k=0...n/2 (rounded down) are retained
    types::Uint nXF = nX, nYF = nY/2 + 1;
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nXF * nYF);
    fftw_plan plan_forward = fftw_plan_dft_r2c_2d(nX, nY, planeDirect.data(), out, FFTW_ESTIMATE);


    fftw_execute(plan_forward); /* repeat as needed */

    Array<types::Complex,2> planeFourier(
        reinterpret_cast<types::Complex *>(out),
        shape(nXF, nYF),
        neverDeleteData);

    t = clock();
    // Calculate the exponential prefactors.
    // Multiplication of Fourier coefficients with this prefactors
    // propagates them to next z plane.i
    // SI units energy: 2 m E/hbar^2 a0 = 2 E/Ha
    
    // the prefectors would only need dimensions (nX/2 +1, nY/2 +1)
    // since they are the same for G and -G. However for the sake of
    // simplicity of calculation, we prepare it here for both.
    // notice that the order of storage is 0...G -G...0
    Array<types::Real,2> prefactors(nXF, nYF);
    types::Real energyTerm = 2 * wfn.energy;
    types::Real deltaZ = wfn.grid.directions[2].incrementVector[2];
    la::Cell cell = la::Cell(wfn.grid.directions);
    vector<types::Real> X = cell.vector(0);
    vector<types::Real> Y = cell.vector(1);
    types::Real dKX = 2 * M_PI / X[0];
    types::Real dKY = 2 * M_PI / Y[1];

    prefactors( Range(0, nXF/2 - 1), Range::all())= exp(- sqrt(tensor::i * dKX * tensor::i * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);
    // tensor::i always starts from 0
    prefactors( Range(nXF/2, nXF - 1), Range::all())= exp(- sqrt( (nXF/2 - tensor::i) * dKX * (nXF/2 - tensor::i) * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);

    
    std::cout << prefactors;
    std::cout << planeDirect;
    std::cout << planeFourier;
    // Sequentially update the cube file
    Array<types::Real,2> tempDirect(nX, nY);
    Array<types::Complex,2> tempFourier(nXF, nYF);
    fftw_plan plan_backward;
    for(int zIndex = startIndex + 1; zIndex <= endIndex; ++zIndex) {
        // Propagate fourier coefficients
        // The y-dim is aready cut in half, for x we have to do it
        planeFourier = prefactors(tensor::i, tensor::j) * planeFourier(tensor::i, tensor::j);
    // The c2r transform destroys its input
        tempFourier = planeFourier;
        plan_backward = fftw_plan_dft_c2r_2d(nX, nY, (fftw_complex*) tempFourier.data(), tempDirect.data(), FFTW_ESTIMATE);
        // Do Fourier-back transform
        fftw_execute(plan_backward); /* repeat as needed */
        tempDirect /= nX * nY;
        std::cout << tempDirect;
        // Copy data
        dataArray(Range::all(), Range::all(), zIndex) = tempDirect(Range::all(), Range::all());

    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(out);

    std::cout << "Time to extrapolate : " << (clock() -t)/1000.0 << " ms\n";
    t = clock();
    types::String outCubeFile = "out.";
    outCubeFile += cubeFile;
    wfn.writeCubeFile(outCubeFile);
    std::cout << "Time to write cube : " << (clock() -t)/1000.0 << " ms\n";

    // Produce z profile
    std::string outZProfile = outCubeFile;
    outZProfile += ".zprofile";
    wfn.writeZProfile(outZProfile, "Z profile of squared extrpolated wave function\n");
    wfn.grid.squareValues();
    outZProfile += ".zprofile.squared";
    wfn.writeZProfile(outZProfile, "Z profile of squared extrpolated wave function\n");


    return true;
}
