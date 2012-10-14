#include "formats/cube.hpp"
#include "formats/gnuplot.hpp"
#include "formats/stm.hpp"
#include "atomistic/fundamental.hpp"
#include "types.hpp"
#include "wrappers/lapack.hpp"
#include "io.hpp"
#include <cmath>

#include <ctime>

#include <boost/format.hpp>
#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi_eol.hpp>
#include <boost/lexical_cast.hpp>

#include <blitz/array.h>

#include <fftw3.h>

namespace formats {
namespace stm {

using namespace types;
using namespace wrappers::lapack;

bool StsCube::initialize(){
    if(levels.size() == 0){
        std::cout 
            << "No states for given energy window. STS will not be performed.";
        return false;
    }

    // Steal from first cube file
    std::list< formats::WfnCube >::const_iterator
        cubeIt = levels.begin();
    
    // Get atoms and proper grid
    this->fileName = cubeIt->fileName;
    this->readCubeFile();
    
    this->title  = "STS data (z = energy)";
    this->description  = str(boost::format(
                "Range [%1.3d V,%1.3d V], delta %1.4d V, sigma %1.4d V")
            % eMin % eMax % deltaE % sigma);
    
    if(modeFlag == CONSTANT_Z){
        Real height = boost::lexical_cast<Real>(modeParameter);
    
        // Get z index of STS plane
        std::vector<types::Real> tempvector;
        tempvector.push_back(0);
        tempvector.push_back(0);
        tempvector.push_back(this->topZCoordinate() + height);
        std::vector<types::Uint> indices;
        this->grid.getNearestIndices(tempvector, indices);
        this-> zIndex = indices[2];
        std::cout << "STS will be performed at z-index " 
            << zIndex << "\n";
        this->title += str(boost::format(" at z-plane %3d") % zIndex);
    }
    else if (modeFlag == PROFILE){
        stm = stm::StmCube();
        std::cout << "Reading z profile from " 
            << modeParameter << "\n";
        stm.readIgorFile(modeParameter);
        this->title += " on z profile read from ";
        this->title += modeParameter;
    }

    // Adjust z dimension for energy
    Uint newIncrementCount =(unsigned int) ((eMax - eMin)/deltaE) + 1;

    if(this->cubeZ == 0){
        // Per default we simply keep the old z extent
        this->grid.directions[2].scaleVector(
                Real(grid.directions[2].getNElements())/
                Real(newIncrementCount));
        this->grid.directions[2] = 
            la::Direction(grid.directions[2].getIncrementVector(),
                    newIncrementCount);
    }
    else{
        // Else the new extent is cubeZ atomic units
        Real newDZ = this->cubeZ / newIncrementCount;
        std::vector<Real> v; v.push_back(0); v.push_back(0); v.push_back(newDZ);
        this->grid.directions[2] = la::Direction(v, newIncrementCount);
        this->grid.originVector[2] = 
            eMin / (eMax - eMin) * cubeZ;
    }
    grid.data = std::vector<Real>(grid.countPoints(), 0.0);

    return true;
}


bool StsCube::calculate(){

    if( ! this->initialize() ) return false;

    std::list< formats::WfnCube >::const_iterator
        cubeIt = levels.begin(),
    cubeEnd = levels.end();

    // Get all necessary planes
    formats::WfnCube tempCube;
    while(cubeIt != cubeEnd){
        std::cout<< "Processing " << cubeIt->fileName << "\n";
        tempCube = *cubeIt;
        tempCube.readCubeFile();
        if(!this->psiSquared) tempCube.squareValues();
        std::vector<Real> tempPlane;

        if(modeFlag == CONSTANT_Z){       
           tempCube.getZPlane(zIndex, tempPlane);
        }
        else if(modeFlag == PROFILE){
            this->interpolateOnZProfile(tempCube.grid, tempPlane);
        }
        
        addLevel(tempPlane, cubeIt->energy); 
        ++cubeIt;
    }
        
    std::cout<< "Done processing cube files...\n\n";

    // Normalize sum, s.th. values are, on average, 1
    //grid *= Real(grid.countPoints()) / grid.sum();

    return true;
}

void StsCube::interpolateOnZProfile(
        const formats::CubeGrid &grid,
        std::vector<Real> &result) const { 
    Uint nX = grid.directions[0].getNElements();
    Uint nY = grid.directions[1].getNElements();
    result.reserve(nX * nY);

    Real dZ = grid.directions[2].getIncrementVector()[2];
    std::vector<Real>::const_iterator zIt = this->stm.stm.begin();
    Real valLow, valHigh, deltaLow;
    Uint zLow;
    for(Uint x = 0; x < nX; ++x){
        for(Uint y = 0; y < nY; ++y){
            zLow = Uint(*zIt / dZ);
            deltaLow = (*zIt - zLow * dZ) / dZ; // Between 0 and 1
            valLow = grid.getDataPoint(x, y, zLow);
            valHigh = grid.getDataPoint(x, y, zLow + 1);

            result.push_back(valLow * (1.0 - deltaLow) + valHigh * deltaLow);

            ++zIt;
        }
    }
}


void StsCube::addLevel(const std::vector<types::Real> &plane,
               types::Real energy){
    
    Uint nEnergies = grid.directions[2].getNElements();

    // Gaussian stuff
    // \$ \frac{1}{\sigma \sqrt{2\pi} e^{-\frac{(x-\mu)^2}{2\sigma^2}} \$
    // = a e^{c(x-b)^2}
    types::Real a = 1.0/(sigma * std::sqrt(2 * M_PI));
    types::Real c = -1.0/(2 * sigma * sigma);
   
     std::vector<Real>::const_iterator 
         planeIt = plane.begin(),
         planeEnd = plane.end();
     std::vector<Real>::iterator 
         dataIt = grid.data.begin();
     while(planeIt != planeEnd){
         for(Uint z = 0; z < nEnergies; ++z){
             Real e = eMin + deltaE * z; 
             *dataIt += *planeIt * a * std::exp(c * (energy - e)*(energy - e));
             ++dataIt;
         }
         ++planeIt;
     }


}

void StmCube::setIsoValue(types::Real isoValue){
    this->setValue = isoValue;
    this->mode = CONSTANT_CURRENT;
    this->stm.clear();
    getZIsoSurface(isoValue, this->stm);
}

void StmCube::setZValue(types::Real zValue){
    this->setValue = zValue;
    this->mode = CONSTANT_HEIGHT;
    this->stm.clear();

    // Find largest z coordinate of atoms
    std::vector<types::Real> tempvector;
    tempvector.push_back(0);
    tempvector.push_back(0);
    tempvector.push_back(this->topZCoordinate() + zValue);
    std::vector<types::Uint> indices;
    this->grid.getNearestIndices(tempvector, indices);
    Uint zIndex = indices[2];
    std::cout << "STM will be performed at z-index " 
        << zIndex << "\n";
    this->getZPlane(zIndex, this->stm);
    std::cout << "STM will be performed at z-index " 
        << zIndex << "\n";
}

bool StmCube::writeIgorFile(String fileName) const {
    Uint nX = grid.directions[0].getNElements();
    Uint nY = grid.directions[1].getNElements();
    std::vector<Real>::const_iterator stmIt = stm.begin();

    Stream result = "";
    for(Uint x = 0; x<nX; ++x){
        for(Uint y=0; y<nY; ++y){
            result += str(boost::format("%12.6e") % *stmIt);
            result += " ";
            ++stmIt;
        }
        result += "\n";
    }

    return io::writeStream(fileName, result);
}


bool StmCube::readIgorFile(String fileName) {
    String content;
    io::readFile(fileName, content);
    
    std::string::const_iterator it = content.begin(),
        end = content.end();

    using boost::spirit::double_;
    using boost::spirit::qi::eol;
    using boost::spirit::qi::parse;
    using boost::spirit::ascii::space;
    if (! parse(
        it,
        end,
        double_ % space,
        stm
        )) throw types::parseError() << types::errinfo_parse("title or description");

   return true;
} 

WfnExtrapolation::WfnExtrapolation(
        std::vector<String> fileNames,
        const cp2k::Spectrum& spectrum,
        Cube&   hartree ,
        Mode          mode    ,
        types::Real   var1    ,
        types::Real   zWidth  ,
        types::Real   approachFrom,
        types::Real   kPlaneMax,
        types::Uint   nLayers
        ){
    
    for(std::vector<String>::const_iterator it = fileNames.begin();
                        it != fileNames.end(); ++it){
        WfnCube tmp;
        tmp.setFileName(*it);
        wfns.push_back(tmp);
    }
        

    this->spectrum = spectrum;
    this->hartree = hartree;
    this->mode    = mode;
    if      ( mode == plane )       this->zStart     = var1;
    else if ( mode == rollingBall ) this->ballRadius = var1;
    else if ( mode == isoSurface )  this->isoValue   = var1;
    this->zWidth  = zWidth;
    this->approachFrom = approachFrom;
    this->kPlaneMax = kPlaneMax;
    this->nLayers     = nLayers;

}

void WfnExtrapolation::execute(){
    using namespace blitz;
    time_t t;

    this->determineRange();
    
    
    // Need some info on the grid
    // SI units energy: 2 m E/hbar^2 a0 = 2 E/Ha
    Uint nX = hartree.nX(), nY = hartree.nY(), nZ = zEndIndex + 1;
    Real dX = hartree.dX(), dY = hartree.dY(), dZ = hartree.dZ();
    
    std::cout << "Original wave vector grid: kx,ky = " 
        << nX << "," << nY << "\n";
    int nKX = nX, nKY = nY, nK = nKX * nKY;
    if (mode == isoSurface || mode == rollingBall){ 
        // May choose to reduce Fourier components
        // limiting oscillation frequency/a.u. to kMax
        Real dKX = 2 * M_PI / (dX * Real(nX));
        Real dKY = 2 * M_PI / (dY * Real(nY));

        Real kXMax = nKX/2.0 * dKX;
        Real kYMax = nKY/2.0 * dKY;
        Real scale = this->kPlaneMax / sqrt(kXMax * kXMax + kYMax * kYMax);
        if( scale < 1.0 ){
            nKX = int( scale * nX);
            nKY = int( scale * nY);
            std::cout << "Reduced wave vector grid: kx,ky = " << nKX << "," << nKY << "\n";
        }

    }
    
    for(std::vector<WfnCube>::iterator wfn = wfns.begin();
            wfn != wfns.end(); ++wfn){

        t = clock();
        std::cout << "--------\n Processing  " << wfn->getFileName() << "\n";
        wfn->readCubeFile(); 
        std::cout << "Time to read cube : " << (clock() -t)/1000.0 << " ms\n";
        
        wfn->grid.resize(nX, nY, nZ);
        Array<types::Real,3> dataArray(
                &(wfn->grid.data[0]),
                shape(nX, nY, nZ),
                neverDeleteData);
        Real E =  spectrum.getLevel(wfn->getSpin(), wfn->getLevel());
        E -= this->surfacePotential;
        wfn->setEnergy(E);
        std::cout << "Energy level " << wfn->getEnergy() << " Ha\n";

        // Produce z profile
        std::string zProfile = io::getFileName(wfn->getFileName());
        WfnCube tmp = WfnCube(*wfn);
        tmp.squareValues();
        zProfile += ".zprofile";
        std::cout << "Writing Z profile to " << zProfile << std::endl;
        tmp.writeZProfile(zProfile, "Z profile of sqared wave function\n");
           

        // Get values on extrapolation surface
        wfn->grid.zSurface(this->zIndices, this->surface);
        types::String s = formats::gnuplot::writeMatrix<Real>(this->surface, nX, nY);
        String sFile = io::getFileName(wfn->getFileName()) + ".exsurf";
        io::writeStream(sFile, s);

        std::cout << "Starting extrapolation ... " << std::endl;
        if ( mode == plane ) 
            this->onPlane(*wfn);
        
        else if ( mode == isoSurface || mode == rollingBall )  
            this->onSurface(*wfn, nKX, nKY);
            
        std::cout << "Time to extrapolate : " << (clock() -t)/1000.0 << " ms\n";

        types::String outCubeFile = "extrapolated.";
        outCubeFile += io::getFileName(wfn->getFileName());
        t = clock();
        std::cout << "Writing extrapolated cube file to " << outCubeFile << std::endl;
        wfn->writeCubeFile(outCubeFile);
        std::cout << "Time to write cube : " << (clock() -t)/1000.0 << " ms\n";

        // Produce z profile
        tmp = *wfn;
        tmp.grid.squareValues();
        std::string outZProfile = outCubeFile;
        outZProfile += ".zprofile";
        std::cout << "Writing Z profile of extrapolated cube file to " << outZProfile << std::endl;
        tmp.writeZProfile(outZProfile, "Z profile of squared extrapolated wave function\n");

        wfn->grid.data.clear();
    }
}

void WfnExtrapolation::determineRange(){
    using namespace blitz;

    Cube tmp(hartree);

    // Some preparations...
    if( mode == plane || mode == rollingBall) {
        // Find highest z-coordinate
        // Note: The slowest index in cube file format is x, so in terms of storage
        //       modifying the number of x-coordinates would be the easiest.
        std::vector< atomistic::Atom >::const_iterator it = hartree.atoms.begin(),
            end = hartree.atoms.end();
        types::Real zTop = it->coordinates[2];
        ++it;
        while(it != end) {
            if( it->coordinates[2] > zTop ) {
                zTop = it->coordinates[2];
            }
            ++it;
        }
    
        // Define region of interpolation
        std::vector<types::Real> tempvector;
        tempvector.push_back(0);
        tempvector.push_back(0);
        tempvector.push_back(zTop + zStart);
        std::vector<types::Uint> indices;
        hartree.grid.getNearestIndices(tempvector, indices);
        this->zStartIndex = indices[2];
        this->zSurfEndIndex = indices[2];
        
        tempvector[2] += zWidth;
        hartree.grid.getNearestIndices(tempvector, indices);
        this->zEndIndex = indices[2];
    }


    // Define the z indices
    if (mode == plane){
        std::cout << "Fourier transform will be performed at z-index " 
            << zStartIndex << "\n";
    
        this->zIndices = std::vector<Uint> (hartree.nX() * hartree.nY(), zStartIndex);
        
    }

    else if ( mode == isoSurface ) {
        if (approachFrom < 0) hartree.grid.zIsoSurfaceOnGrid(isoValue, this->zIndices);
        else                  hartree.grid.zIsoSurfaceOnGrid(isoValue, this->zIndices, approachFrom);
        
        Array<types::Uint,2> zIndicesBlitz(
                &zIndices[0],
                shape(hartree.nX(), hartree.nY()),
                neverDeleteData);
        
        // Write some info about plane 
        zStartIndex = min(zIndicesBlitz);
        Real dZ = hartree.grid.directions[2].getIncrementVector()[2];
        zSurfEndIndex = max(zIndicesBlitz);
        zEndIndex = zSurfEndIndex + Uint(zWidth / dZ);

        std::cout << "Isosurface ranges over " << (zSurfEndIndex-zStartIndex)*dZ
                  << " from z = " << zStartIndex * dZ 
                  << " to " << max(zIndicesBlitz) * dZ << " [a.u.].\n";
        std::cout << "Extrapolation will be extended to "
                  << zEndIndex * dZ << " [a.u.]\n";
    }

    else if (mode == rollingBall){
        Uint nX = hartree.nX(), nY = hartree.nY(), nZ = zEndIndex + 1;
        Real dX = hartree.dX(), dY = hartree.dY(), dZ = hartree.dZ();

        // TODO if (approachFrom > 0) zMax  = approachFrom;


        std::vector<bool> boolGrid(nX*nY*nZ, false);
        std::vector<atomistic::Atom>::const_iterator 
            it  = hartree.atoms.begin(),
            end = hartree.atoms.end();
        int rI = int(ballRadius / dX) + 1;
        int rJ = int(ballRadius / dY) + 1;
        
        Real x,y,z;
        int iX, iY, iZ;
        Real r2 = ballRadius * ballRadius;

        // Build the boolean Grid
        std::vector<types::Uint> indices;
        std::vector<Real> coords;
        while(it != end){
            indices.clear();
            coords = it->getCoordinates();
            hartree.grid.getNearestIndices(coords, indices);

            for(int i = -rI; i <= rI; ++i)
                for(int j = -rJ; j <= rJ; ++j){
                    x = i * dX;
                    y = j * dY;
                    z = r2 - x*x - y*y;

                    if (z > 0){
                        iX = (indices[0] + i) % nX;
                        iY = (indices[1] + j) % nY;
                        iZ = int(indices[2]) + int(sqrt(z) / dZ);
                        if (iZ > zEndIndex) iZ = zEndIndex;

                        boolGrid[iX*nY*nZ + iY*nZ + iZ] = true;
                    }
                }

            ++it;
        }

               

        // Find the topmost 'true' surface
        zIndices = std::vector<Uint>(nX*nY);
        for(int i = 0; i < nX; ++i)
            for(int j = 0; j < nY; ++j){
                int k = nZ;
                while(k > 0){
                    --k;
                    if (boolGrid[i*nY*nZ + j*nZ + k] == true){
                        zIndices[i*nY + j] = k;
                        break;
                    }
                }
                if (k == 0){
                    std::cout << "Warning: Ball rolled to bottom of cube file without finding an atom.\n";
                    zIndices[i*nY + j] = 0;
                }
            }

        
//        Old, much too slow way of doing things
//        std::vector<Real> tmp;
//        for(Uint x = 0; x < nX; ++x){
//            std::cout << "Here\n";
//            for(Uint y = 0; y < nY; ++y)
//                for(Uint z = zMax; z > 0; --z){
//                    tmp.push_back(x*hartree.dX());
//                    tmp.push_back(y*hartree.dY());
//                    tmp.push_back(z*hartree.dZ());
//                    if(hartree.distance(tmp) < ballRadius){
//                        zIndices[y*nX+x] = z-1;
//                        break;
//                    }
//                }
//        }
        
        Array<types::Uint,2> zIndicesBlitz(
                &zIndices[0],
                shape(hartree.nX(), hartree.nY()),
                neverDeleteData);
        
        // Write some info about plane 
        zStartIndex = min(zIndicesBlitz);
        zSurfEndIndex = max(zIndicesBlitz);
        //Real dZ = hartree.dZ();
        zEndIndex = zSurfEndIndex + Uint(zWidth / dZ);

        std::cout << "Rolling-ball surface ranges over " << (zSurfEndIndex-zStartIndex)*dZ
                  << " from z = " << zStartIndex * dZ 
                  << " to " << zSurfEndIndex * dZ << " [a.u.].\n";
        std::cout << "Extrapolation will be extended to "
                  << zEndIndex * dZ << " [a.u.]\n";

    }
   
    // Write extrapolation z indices 
    String s = formats::gnuplot::writeMatrix<Uint>(this->zIndices, hartree.nX(), hartree.nY());
    String sFile = io::getFileName(hartree.getFileName()) + ".zindices";
    io::writeStream(sFile, s);

    this->setSurfacePotential();
    

}
   
void WfnExtrapolation::setSurfacePotential(){
    using namespace blitz;
   

    std::vector<Real> hartreeVector;
    hartree.grid.zSurface(zIndices, hartreeVector);
    Array<Real, 2> hartreeSurface (
            &hartreeVector[0],
            shape(hartree.nX(), hartree.nY()),
            neverDeleteData);
    this->surfacePotential = mean(hartreeSurface);
    Array<Real, 3> hartreeGrid (
            &(hartree.grid.data[0]),
            shape(hartree.nX(), hartree.nY(), hartree.nZ()),
            neverDeleteData);
    Real vacuumPotential = max(hartreeGrid);

    std::cout << "Average Hartree potential on extrapolation surface is " 
              << surfacePotential << " Ha\n";
    std::cout << "This is " << vacuumPotential - surfacePotential 
              << " Ha below the maximum Hartree potential.\n";

    // Print hartree potential for gnuplot
    types::String s = formats::gnuplot::writeMatrix(hartreeVector, hartree.nX(), hartree.nY());
    String surfaceFile = io::getFileName(hartree.getFileName()) + ".exsurf";
    io::writeStream(surfaceFile, s);

}


void WfnExtrapolation::onPlane(WfnCube& wfn){
    using namespace blitz;
    time_t t;
        
    // Need some info on the grid
    // SI units energy: 2 m E/hbar^2 a0 = 2 E/Ha
    Uint nX = wfn.nX(), nY = wfn.nY(), nZ = wfn.nZ();
    Real dX = hartree.dX(), dY = hartree.dY(), dZ = hartree.dZ();
    Real dKX = 2 * M_PI * dZ / (dX * Real(nX));
    Real dKY = 2 * M_PI * dZ / (dY * Real(nY));
    
    Real E = wfn.getEnergy() * dZ * dZ;
    Array<types::Real,3> dataArray(
            &(wfn.grid.data[0]),
            shape(nX, nY, nZ),
            neverDeleteData);

    // Do a real 2 complex fft
    // Since a(-k)=a(k)^* (or in terms of indices: a[n-k]=a[k]^*)
    // only a[k], k=0...n/2+1 (division rounded down) are retained
    types::Uint nXF = nX, nYF = nY/2 + 1;
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nXF * nYF);

    Array<types::Real,2> planeDirect(
            &surface[0], 
            shape(nX, nY), 
            neverDeleteData);

    fftw_plan plan_forward = fftw_plan_dft_r2c_2d(nX, nY, planeDirect.data(), out, FFTW_ESTIMATE);
    fftw_execute(plan_forward);


    Array<types::Complex,2> planeFourier(
            reinterpret_cast<types::Complex *>(out),
            shape(nXF, nYF),
            neverDeleteData);
    
    std::vector<Real> fourierReal(nXF*nYF);
    for(int i =0;i<nXF;++i)
        for(int j =0;j<nYF;++j)
            fourierReal[i*nYF+j] = abs(planeFourier(i,j));

    // Plotting the coefficients of the Fourier
    String s = formats::gnuplot::writeMatrix<Real>(fourierReal, nXF, nYF);
    String hartleyFile = io::getFileName(wfn.getFileName()) + ".fourier";
    io::writeStream(hartleyFile, s);


    // Calculate the exponential prefactors.
    // Multiplication of Fourier coefficients with these prefactors
    // propagates them to next z plane.
    // The prefectors would only need dimensions (nX/2 +1, nY/2 +1)
    // since they are the same for G and -G. However for the sake of
    // simplicity of calculation, we prepare them here for (nX, nY/2+1) = (nXF,nYF).
    t = clock();
    Array<types::Real,2> prefactors(nXF, nYF);

    // Notice that the order of storage is 0...G -G...-1 for uneven nX
    // and 0...G-1 G -(G-1)...-1 for even NX
    prefactors( Range(0, nXF/2), Range::all())= exp(- sqrt(tensor::i * dKX * tensor::i * dKX + tensor::j *dKY * tensor::j * dKY - 2 * E) );
    // tensor::i always starts from 0, i.e. it ranges from 0...nXF/2-1 (nX even) or 0...nXF/2 (nX uneven)
    if(nX % 2 == 1)
        prefactors( Range(nXF/2 + 1, nXF - 1), Range::all())= exp(- sqrt( (nXF/2 - tensor::i) * dKX * (nXF/2 - tensor::i) * dKX + tensor::j *dKY * tensor::j * dKY - 2 * E) );
    else
        prefactors( Range(nXF/2 + 1, nXF - 1), Range::all())= exp(- sqrt( (nXF/2 -1 - tensor::i) * dKX * (nXF/2 - 1 - tensor::i) * dKX + tensor::j *dKY * tensor::j * dKY - 2 * E) );

    // Sequentially update the cube file
    Array<types::Real,2> tempDirect(nX, nY);
    Array<types::Complex,2> tempFourier(nXF, nYF);
    fftw_plan plan_backward;
    for(int zIndex = zStartIndex + 1; zIndex < nZ; ++zIndex) {
        // Propagate fourier coefficients
        planeFourier = prefactors(tensor::i, tensor::j) * planeFourier(tensor::i, tensor::j);

        // The c2r transform destroys its input
        tempFourier = planeFourier;
        plan_backward = fftw_plan_dft_c2r_2d(nX, nY, (fftw_complex*) tempFourier.data(), tempDirect.data(), FFTW_ESTIMATE);


        // Do Fourier-back transform
        fftw_execute(plan_backward); /* repeat as needed */
        tempDirect /= nX * nY;

        // Copy data
        dataArray(Range::all(), Range::all(), zIndex) = tempDirect(Range::all(), Range::all());
    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(out);
}

void WfnExtrapolation::onSurface(WfnCube& wfn, Uint nKX, Uint nKY){
    using namespace blitz;
    time_t t;
    
    // Need some info on the grid
    // SI units energy: 2 m E/hbar^2 a0 = 2 E/Ha
    Uint nX = hartree.nX(), nY = hartree.nY(), nZ = zEndIndex + 1;
    Real dX = hartree.dX(), dY = hartree.dY(), dZ = hartree.dZ();
    Real dKX = 2 * M_PI * dZ / (dX * Real(nX));
    Real dKY = 2 * M_PI * dZ / (dY * Real(nY));
    
    Array<types::Real,3> dataArray(
            &(wfn.grid.data[0]),
            shape(nX, nY, nZ),
            neverDeleteData);
    Real E = wfn.getEnergy() * dZ * dZ;

    int nK = nKX * nKY;
   
    // Instead of a wave-number based criterion that is the
    // same for all wave functions, we may also choose some adaptive
    // criterion such as the one below (upper limit for z-decay of
    // involved basis functions.
//    // May choose to reduce k-grid 
//    Real kXMax = nKX/2.0 * dKX;
//    Real kYMax = nKY/2.0 * dKY;
//    // Take only low-freq. Fourier components such that no basis function
//    // decays faster than f(z) \propto 10^{-minExpGoal * z}
//    Uint zSpan = zSurfEndIndex - zStartIndex;
//
//    // Translate decayCutoff to minimum exponent per z-index
//    Real minExpGoal = this->decayCutoff * dZ * log(10.0);   
//    Real scaling = sqrt( (minExpGoal * minExpGoal + 2*E)
//            / (kXMax * kXMax + kYMax * kYMax) );
//    if(scaling < 1.0){
//        nKX = int( scaling * nX); kXMax = nKX/2.0 * dKX;
//        nKY = int( scaling * nY); kYMax = nKY/2.0 * dKY;
//        std::cout << "Reduced wave vector grid: kx,ky = " << nKX << "," << nKY << "\n";
//    }
//    Real minExp = - sqrt( -2.0 * E + kXMax * kXMax + kYMax * kYMax);
//    std::cout << "Strongest considered decay: x " << exp(minExp/dZ) << " per a.u.\n";
//    minExp *= Real(zSurfEndIndex) - Real(zStartIndex);
//    std::cout << "Smallest exponential term in matrix: " << exp(minExp) << "\n";
        
    

    // Build matrix for modified discrete Hartley transform
    Uint n = nX * nY;
    std::cout << "Matrix dimension is " << n << "x" << nK 
        << ", i.e. " << Real(nLayers*n*nK*8)  / (1024.0 * 1024.0) << " MBytes per matrix\n";
    //exit(0); 
    //        // Test: Extrapolation on plane
    //        for(int i = 0; i < zIndices.size();++i){
    //            zIndices[i] = 128;
    //        }
    //        this->zSurfEndIndex = 128;
    //        this->zStartIndex = 128;
    //        E = -0.160023 * dZ * dZ;

    std::vector<Real> A(n*nK*nLayers);
    int iX, iY, jX, jY;
    Real kX, kY, kZ;
    Real arg;
    t = clock();


    for(Uint j = 0; j < nK; ++j){
        jX = j / nKY; jX = (2* jX < nKX) ? jX : jX - nKX;
        jY = j % nKY; jY = (2* jY < nKY) ? jY : jY - nKY;

        kX = jX * dKX;
        kY = jY * dKY;
        kZ = sqrt( -2.0 * E + kX * kX+ kY * kY);

        for(Uint i = 0; i < n; ++i){
            iX = i / nY; iY = i % nY;
            arg = 2.0 * M_PI * ( Real(iX * jX) / Real(nX)  + Real(iY * jY) / Real(nY) );

            for(Uint l = 0; l < nLayers; ++l){
                // Need transposed matrix for Fortran
                A[l*n+i + j*n*nLayers] = (sin(arg) + cos(arg)) *
                    exp(- kZ * (int(zIndices[i]) - int(l) - int(zStartIndex) ) );
            }
        }
    }
    std::cout << "Matrix created in : " << (clock() -t)/1000.0 << " ms\n";


    // Prepare wave function values
    int iZ;
    std::vector<Real> hartley(n*nLayers);
    for(Uint l = 0; l < nLayers; ++l){
        for(Uint i = 0; i < n; ++i){
            iX = i/nY;
            iY = i%nY;
            iZ = int(zIndices[i]) - int(l);
            hartley[l*n + i] = wfn.grid.data[iX*nY*nZ+ iY*nZ + iZ];
        }
    }


    // Solve the linear system
    char TRANS_ = 'N';
    int  M_     = n*nLayers;
    int  N_     = nK;
    int  NRHS_  = 1;
    std::vector<Real> S(N_);
    Real RCOND_ = 1e-5;
    int  RANK_;
    std::vector<Real> WORK(1);
    int  LWORK_ = -1;
    std::vector<int> IWORK(1);
    int  ILWORK_ = -1;
    int  INFO_  = 0;
    // Writes correct WORK dimensions to WORK[0], IWORK[0]
    // SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
    //                   WORK, LWORK, IWORK, INFO )
    dgelsd_(&M_, &N_, &NRHS_, 
            &*A.begin(), &M_, 
            &*hartley.begin(), &M_, 
            &*S.begin(), &RCOND_, &RANK_,
            &*WORK.begin(), &LWORK_, &*IWORK.begin(),
            &INFO_);
    LWORK_ = int(abs(WORK[0])); 
    WORK =  std::vector<Real>(LWORK_);
    ILWORK_ = int(abs(IWORK[0])); 
    IWORK =  std::vector<int>(ILWORK_);

    // Now solve this little puzzle
    std::cout << "Starting optimization of overdetermined system\n"; 
    t = clock();

    dgelsd_(&M_, &N_, &NRHS_, 
            &*A.begin(), &M_, 
            &*hartley.begin(), &M_, 
            &*S.begin(), &RCOND_, &RANK_,
            &*WORK.begin(), &LWORK_, &*IWORK.begin(),
            &INFO_);
    //        dgels_(&TRANS_, &M_, &N_, &NRHS_, 
    //                &*A.begin(), &M_, 
    //                &*hartley.begin(), &M_, 
    //                &*WORK.begin(), &LWORK_,
    //                &INFO_);
    if( INFO_ != 0 ){
        throw types::runtimeError() 
            << types::errinfo_runtime("LAPACK Error: Problem with solution of least-squares problem.");
    }
    std::cout << "Found solution in : " << (clock() -t)/1000.0 << " ms\n";

    // Print condition number of problem
    std::cout << "Condition number: " << S[0] / S[nK -1] << "\n";

    //        // Plotting the coefficients of the hartley transform
    //        types::String s = formats::gnuplot::writeMatrix<Real>(hartley, nKX, nKY);
    //        String hartleyFile = io::getFileName(wfn->getFileName()) + ".hartley";
    //        io::writeStream(hartleyFile, s);

    // Re-layout hartley components to nX-nY
    std::vector<Real> tmp = hartley;
    hartley = std::vector<Real>(nX*nY, 0.0);
    Real i2, j2;
    for(int i = 0; i < nKX; ++i){
        i2 = (2 * i < nKX) ? i :  int(nX) + i - nKX ;
        for(int j = 0; j < nKY; ++j){
            j2 = (2 * j < nKY) ? j : int(nY) + j - nKY;
            hartley[i2 * nY + j2] = tmp[i* nKY + j];
        }
    }

    //        // Plotting the coefficients of the reformed hartley transform
    //        s = formats::gnuplot::writeMatrix<Real>(hartley, nX, nY);
    //        hartleyFile = io::getFileName(wfn->getFileName()) + ".hartley2";
    //        io::writeStream(hartleyFile, s);

    // Producing Fourier transform
    int nXF = nX; int nYF = nY/2 +1;
    std::vector<Complex> fourier(nXF*nYF);
    int iT,jT;
    for(int i = 0; i < nXF; ++i){
        for(int j = 0; j < nYF; ++j){
            iT = (int(nX) - i) % nX;
            jT = (int(nY) - j) % nY;

            fourier[i*nYF+j] = 
                1/2.0 * (               hartley[i*nY+j] + hartley[iT*nY+jT]  
                        - Complex(0,1.0)*( hartley[i*nY+j] - hartley[iT*nY+jT]) );
        }
    }

    
    // Plotting the coefficients of the Fourier transform
    std::vector<Real> fourierReal(fourier.size());
    for(int i =0;i<fourierReal.size();++i)
        fourierReal[i] = abs(fourier[i]);
    String s = formats::gnuplot::writeMatrix<Real>(fourierReal, nXF, nYF);
    String fileName = io::getFileName(wfn.getFileName()) + ".fourier";
    io::writeStream(fileName, s);


    // Calculating the exponential prefactors to propagate coefficients to
    // next z plane
    t = clock();
    std::vector<Real> pref(nXF*nYF);    
    for(int i = 0; i < nXF; ++i){
        kX = (2*i < nX) ? i : i - int(nX) ;
        kX *= dKX;
        for(int j = 0; j < nYF; ++j){
            kY = (2*j < nY) ? j : j - int(nY);
            kY *= dKY;
            kZ = sqrt( -2.0 * E + kX * kX+ kY * kY);
            pref[i*nYF + j] = exp(- kZ );
        }
    }


    Array<types::Real,2> prefactors(
            &pref[0],
            shape(nXF, nYF),
            neverDeleteData);

    Array<types::Complex,2> planeFourier(
            &fourier[0],
            shape(nXF, nYF),
            neverDeleteData);

    Array<types::Complex,2> tempFourier(nXF, nYF);
    Array<types::Real,2> tempDirect(nX, nY);
    fftw_plan plan_backward;

    for(int zIndex = zStartIndex+1; zIndex < nZ; ++zIndex) {
        // Propagate hartley coefficients (could use matmul here)
        planeFourier = prefactors(tensor::i, tensor::j) * planeFourier(tensor::i, tensor::j);

        // The c2r transform destroys its input
        tempFourier = planeFourier;
        plan_backward = fftw_plan_dft_c2r_2d(
                nX, nY, 
                (fftw_complex*) tempFourier.data(), 
                tempDirect.data(), 
                FFTW_ESTIMATE);

        // Perform Fourier transform backwards
        fftw_execute(plan_backward);

        // Copy data
        if (zIndex > zSurfEndIndex) 
            dataArray(Range::all(), Range::all(), zIndex) = tempDirect(Range::all(), Range::all());
        else {
            for(int x = 0; x < nX; ++x){
                for(int y = 0; y < nY; ++y){
                    if (zIndex > zIndices[x*nY + y]) {
                        dataArray(x, y, zIndex) = tempDirect(x, y);
                    }
                }
            }
        }

    }
    fftw_destroy_plan(plan_backward);

}

}
}
