#include "formats/cube.hpp"
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

void StmCube::setIsoLevel(types::Real isoValue){
    this->stm.clear();
    getZIsoSurface(isoValue, this->stm);
    this->isoLevel = isoValue;
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
        types::String fileName,
        const cp2k::Spectrum& spectrum,
        const Cube&   hartree ,
        Mode          mode    ,
        types::Real   var1    ,
        types::Real   zWidth  ){
    wfn     = WfnCube();  
    time_t t = clock();
    wfn.readCubeFile(fileName.c_str()); 
    std::cout << "Time to read cube : " << (clock() -t)/1000.0 << " ms\n";
    wfn.setEnergy( spectrum.getLevel(wfn.getSpin(), wfn.getLevel()) );

    this->hartree = &hartree;
    this->mode    = mode;
    if ( mode == constantZ ) this->zStart   = var1;
    else                     this->isoValue = var1;
    this->zWidth  = zWidth;
}

void WfnExtrapolation::execute(){
    using namespace blitz;
    time_t t;

    std::cout << "--------\n Extrapolating " << wfn.getFileName() << "\n";
 
    this->determineRange();
    this->adjustEnergy();
    
    types::Uint nX = wfn.nX(), nY = wfn.nY(), nZ = zEndIndex;
    wfn.grid.resize(nX, nY, nZ);
    Array<types::Real,3> dataArray(
        &wfn.grid.data[0],
        shape(nX, nY, nZ),
        neverDeleteData);


    if ( mode == constantZ ){
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

        // Calculate the exponential prefactors.
        // Multiplication of Fourier coefficients with these prefactors
        // propagates them to next z plane.
        // The prefectors would only need dimensions (nX/2 +1, nY/2 +1)
        // since they are the same for G and -G. However for the sake of
        // simplicity of calculation, we prepare them here for (nX, nY/2+1) = (nXF,nYF).
        t = clock();
        Array<types::Real,2> prefactors(nXF, nYF);
        // SI units energy: 2 m E/hbar^2 a0 = 2 E/Ha
        types::Real energyTerm = 2 * wfn.energy;
        std::cout << "Energy level " << wfn.energy << " Ha\n";
        types::Real deltaZ = wfn.grid.directions[2].getIncrementVector()[2];
        la::Cell cell = la::Cell(wfn.grid.directions);
        vector<types::Real> X = cell.vector(0);
        vector<types::Real> Y = cell.vector(1);
        types::Real dKX = 2 * M_PI / X[0];
        types::Real dKY = 2 * M_PI / Y[1];

        // Notice that the order of storage is 0...G -G...-1 for uneven nX
        // and 0...G-1 G -(G-1)...-1 for even NX
        prefactors( Range(0, nXF/2), Range::all())= exp(- sqrt(tensor::i * dKX * tensor::i * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);
        // tensor::i always starts from 0, i.e. it ranges from 0...nXF/2-1 (nX even) or 0...nXF/2 (nX uneven)
        if(nX % 2 == 1)
            prefactors( Range(nXF/2 + 1, nXF - 1), Range::all())= exp(- sqrt( (nXF/2 - tensor::i) * dKX * (nXF/2 - tensor::i) * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);
        else
            prefactors( Range(nXF/2 + 1, nXF - 1), Range::all())= exp(- sqrt( (nXF/2 -1 - tensor::i) * dKX * (nXF/2 - 1 - tensor::i) * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);

        // Sequentially update the cube file
        Array<types::Real,2> tempDirect(nX, nY);
        Array<types::Complex,2> tempFourier(nXF, nYF);
        fftw_plan plan_backward;
        for(int zIndex = zStartIndex + 1; zIndex <= zEndIndex; ++zIndex) {
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
    else if ( mode == isoSurface ) {
        
        // Build matrix for modified Fourier transform
        Uint n = nX * nY;
        Uint iX, iY;
        Real kX, kY, kZ;

        std::vector<Complex> A;
        for(Uint i = 0; i < n; ++i){
            for(Uint j = 0; j < n; ++j){

                iX = i & nY; iY = i / nY;
                if( 2* iX > nX) kX = (nX - i)/ Real( nX );
                else kX = i / Real( nX );
                if( 2* iY > nY) kY = (nY - i)/ Real( nY );
                else kY = i / Real( nY );

                kZ = sqrt( 2 * wfn.getEnergy() + kX * kX + kY * kY);
               
                // Need transposed matrix for Fortran lapack
                A[i + j*n] = exp( 
                            2.0 * M_PI * Complex(0,1) * 
                            ( Real(i % nY) * Real(j % nY) + 
                              Real(i / nY) * Real(j / nY) ) / Real(n)
                            - kZ * (zIndices[i] - zStartIndex) );

            }

        }
        
        // Solve the linear system, overwriting values with Fourier coefficients
        int nVectors = 1; int info; int N = n;
        std::vector<int> permutations(n);  
        
        std::vector<Complex> fourierCoefficients(n);
        for(Uint i = 0; i < surface.size(); ++i) fourierCoefficients[i] = surface[i];
        std::cout << "Solving linear system of dimension " << n 
                  << ", i.e. " << n*n*16 /1024.0/ 1024.0 << " GBytes per matrix\n";
        t = clock();
        zgesv_(&N, &nVectors, &*A.begin(), &N, &*permutations.begin(), &*fourierCoefficients.begin(), &N, &info);
        if( info != 0 ) throw types::runtimeError() 
                           << types::errinfo_runtime("Linear system could not be solved.");
       
        // Now, the resulting coefficients should hopefully be hermitian. We have to check.
        Real normTot, normNonH; 
        Complex aIJ, aJI;
        for(Uint x = 0; x < nX; ++x){
            for(Uint y = 0; y < nY; ++y){
                aIJ = fourierCoefficients[x*nY + y];
                aJI = fourierCoefficients[y*nY + x];
                 
                normTot += abs(aIJ *aIJ);
                normNonH += abs( (aIJ - conj(aJI)) * (aIJ - conj(aJI)));
            }
        }

        std::cout << "Frobenius norm of Fourier coefficients: " << sqrt(abs(normTot)) <<
                     ", Forbenis norm of A-A^H: " << sqrt(abs(normNonH)) << "\n";


//        // Calculate the prefactors
//        // Calculate the exponential prefactors.
//        // Multiplication of Fourier coefficients with these prefactors
//        // propagates them to next z plane.
//        // The prefectors would only need dimensions (nX/2 +1, nY/2 +1)
//        // since they are the same for G and -G. However for the sake of
//        // simplicity of calculation, we prepare them here for (nX, nY/2+1) = (nXF,nYF).
//        t = clock();
//        Array<types::Real,2> prefactors(nXF, nYF);
//        // SI units energy: 2 m E/hbar^2 a0 = 2 E/Ha
//        types::Real energyTerm = 2 * wfn.energy;
//        std::cout << "Energy level " << wfn.energy << " Ha\n";
//        types::Real deltaZ = wfn.grid.directions[2].getIncrementVector()[2];
//        la::Cell cell = la::Cell(wfn.grid.directions);
//        vector<types::Real> X = cell.vector(0);
//        vector<types::Real> Y = cell.vector(1);
//        types::Real dKX = 2 * M_PI / X[0];
//        types::Real dKY = 2 * M_PI / Y[1];
//
//        // Notice that the order of storage is 0...G -G...-1 for uneven nX
//        // and 0...G-1 G -(G-1)...-1 for even NX
//        prefactors( Range(0, nXF/2), Range::all())= exp(- sqrt(tensor::i * dKX * tensor::i * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);
//        // tensor::i always starts from 0, i.e. it ranges from 0...nXF/2-1 (nX even) or 0...nXF/2 (nX uneven)
//        if(nX % 2 == 1)
//            prefactors( Range(nXF/2 + 1, nXF - 1), Range::all())= exp(- sqrt( (nXF/2 - tensor::i) * dKX * (nXF/2 - tensor::i) * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);
//        else
//            prefactors( Range(nXF/2 + 1, nXF - 1), Range::all())= exp(- sqrt( (nXF/2 -1 - tensor::i) * dKX * (nXF/2 - 1 - tensor::i) * dKX + tensor::j *dKY * tensor::j * dKY - energyTerm) * deltaZ);
//
//        fftw_plan plan_backward;
//        for(int zIndex = zStartIndex + 1; zIndex <= zEndIndex; ++zIndex) {
//            // Propagate fourier coefficients
//            planeFourier = prefactors(tensor::i, tensor::j) * planeFourier(tensor::i, tensor::j);
//
//            // The c2r transform destroys its input
//            tempFourier = planeFourier;
//            plan_backward = fftw_plan_dft_c2r_2d(nX, nY, (fftw_complex*) tempFourier.data(), tempDirect.data(), FFTW_ESTIMATE);
//
//            // Do Fourier-back transform
//            fftw_execute(plan_backward); /* repeat as needed */
//            tempDirect /= nX * nY;
//
//            // Copy data
//            dataArray(Range::all(), Range::all(), zIndex) = tempDirect(Range::all(), Range::all());
//        }
//
//        fftw_destroy_plan(plan_forward);
//        fftw_destroy_plan(plan_backward);
    }



    //    // Produce Fourier coefficient for gnuplot    
    //    types::Complex *pt =  reinterpret_cast<types::Complex *>(out);
    //    for(types::Uint i = 0; i < nXF; ++i){
    //            for(types::Uint j = 0;j< nYF; ++j){
    //                    std::cout << std::abs(*pt) << " ";
    //                    ++pt;
    //            }
    //            std::cout << std::endl;
    //    }


    std::cout << "Time to extrapolate : " << (clock() -t)/1000.0 << " ms\n";

}

void WfnExtrapolation::determineRange(){
    using namespace blitz;

    WfnCube tmp = wfn;

    if( mode == constantZ ) {
        // Find highest z-coordinate
        // Note: The slowest index in cube file format is x, so in terms of storage
        //       modifying the number of x-coordinates would be the easiest.
        std::vector< atomistic::Atom >::const_iterator it = wfn.atoms.begin(),
            end = wfn.atoms.end();
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
        wfn.grid.getNearestIndices(tempvector, indices);
        this->zStartIndex = indices[2];
        std::cout << "Fourier transform will be performed at z-index " 
            << zStartIndex << "\n";
        tempvector[2] += zWidth;
        wfn.grid.getNearestIndices(tempvector, indices);
        this->zEndIndex = indices[2];
    
        // Get z-plane for interpolation
        Array<types::Real,3> dataArray(
                &wfn.grid.data[0],
                shape(wfn.nX(), wfn.nY(), wfn.nZ()),
                neverDeleteData);
        Array<types::Real,2> plane = dataArray(Range::all(), Range::all(), zStartIndex);
        this->surface.assign(plane.begin(), plane.end());
    }

    else if ( mode == isoSurface ) {
        tmp.abs();
        tmp.grid.zIsoSurfaceOnGrid(isoValue, this->zIndices, this->surface);
        Array<types::Uint,2> zIndicesBlitz(
                &zIndices[0],
                shape(wfn.nX(), wfn.nY()),
                neverDeleteData);
    
    
        zStartIndex = min(zIndicesBlitz);
        Real dZ = wfn.grid.directions[2].getIncrementVector()[2];
        zEndIndex = max(zIndicesBlitz) + Uint(zWidth / dZ);

        std::cout << "Isosurface ranges from z = " << zStartIndex * dZ 
                  << " to " << max(zIndicesBlitz) * dZ << " [a.u.]\n";
        std::cout << "Interpolation will be extended to "
                  << zEndIndex * dZ << " [a.u.]\n";
    }
    
    // Produce z profile
    std::string zProfile = io::getFileName(wfn.getFileName());
    tmp.squareValues();
    zProfile += ".zprofile";
    std::cout << "Writing Z profile to " << zProfile << std::endl;
    tmp.writeZProfile(zProfile, "Z profile of sqared wave function\n");

}
   
void WfnExtrapolation::adjustEnergy(){
    using namespace blitz;

    Real vacuumPotential = 0;

    if (mode == constantZ ) {
        // Get z-plane for interpolation
        Array<types::Real,3> dataArray(
                const_cast<Real*>(&(hartree->grid.data[0])),
                shape(hartree->nX(), hartree->nY(), hartree->nZ()),
                neverDeleteData);
        Array<types::Real,2> hartreePlane = dataArray(Range::all(), Range::all(), zStartIndex);
        
        // Shift levels by "vacuum" Hartree potential
        // (aveaged value at zStartIndex)
        vacuumPotential = mean(hartreePlane);
        std::cout << "Vacuum hartree potential is " << vacuumPotential << " Ha\n";
        
        hartreePlane    = dataArray(Range::all(), Range::all(), zEndIndex);
        Real hartreeEnd = mean(hartreePlane);
        std::cout << "Hartree potential changes by " << 
            hartreeEnd - vacuumPotential << " Ha over extrapolation region.\n";

//        // Write z plane of hartree potential for gnuplot
//        std::vector<types::Real> hartreePlane;
//        hartree.getZPlane(zStartIndex, hartreePlane);
//        types::String s = formats::gnuplot::writeMatrix(
//                hartreePlane,
//                hartree.grid.directions[0].getNElements(),
//                hartree.grid.directions[1].getNElements());
//        io::writeStream("hartree.zplane",s);
    }
    else if (mode == isoSurface){
        std::vector<Real> hartreeValues;
        const_cast<Cube *>(hartree)->grid.zSurface(zIndices, hartreeValues);
        Array<types::Real,2> hartreeBlitz(
                &hartreeValues[0],
                shape(hartree->nX(), hartree->nY()),
                neverDeleteData);
        vacuumPotential = mean(hartreeBlitz);
        std::cout << "Vacuum hartree potential is " << vacuumPotential << " Ha\n";
    }
    
    wfn.setEnergy( wfn.getEnergy()- vacuumPotential);

}

void WfnExtrapolation::writeWfnCube() const {
    types::String outCubeFile = "extrapolated.";
    outCubeFile += io::getFileName(wfn.getFileName());
    writeWfnCube(outCubeFile);
}
void WfnExtrapolation::writeWfnCube(types::String fileName) const {
    time_t t = clock();
    std::cout << "Writing extrapolated cube file to " 
        << fileName << std::endl;
    wfn.writeCubeFile(fileName);
    std::cout << "Time to write cube : " << (clock() -t)/1000.0 << " ms\n";

    // Produce z profile
    WfnCube temp = wfn;
    temp.grid.squareValues();
    std::string outZProfile = fileName;
    outZProfile += ".zprofile";
    std::cout << "Writing Z profile of extrapolated cube file to " 
        << outZProfile << std::endl;
    temp.writeZProfile(outZProfile, "Z profile of squared extrpolated wave function\n");
}


}
}
