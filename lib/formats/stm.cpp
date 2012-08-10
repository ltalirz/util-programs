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
        types::Real   zWidth  ,
        types::Real   approachFrom,
        types::Real   decayCutoff,
        types::Uint   nLayers
        ){
    wfn     = WfnCube();  
    time_t t = clock();
    std::cout << "--------\n Processing  " << fileName << "\n";
    wfn.readCubeFile(fileName.c_str()); 
    std::cout << "Time to read cube : " << (clock() -t)/1000.0 << " ms\n";
    wfn.setEnergy( spectrum.getLevel(wfn.getSpin(), wfn.getLevel()) );

    this->hartree = &hartree;
    this->mode    = mode;
    if ( mode == constantZ ) this->zStart   = var1;
    else                     this->isoValue = var1;
    this->zWidth  = zWidth;
    this->approachFrom = approachFrom;
    this->decayCutoff = decayCutoff;
    this->nLayers     = nLayers;

}

void WfnExtrapolation::execute(){
    using namespace blitz;
    time_t t;

    this->determineRange();
    this->adjustEnergy();
    
    types::Uint nX = wfn.nX(), nY = wfn.nY(), nZ = zEndIndex + 1;
    wfn.grid.resize(nX, nY, nZ);
    
    // Need some info on the grid
    // SI units energy: 2 m E/hbar^2 a0 = 2 E/Ha
    la::Cell cell = la::Cell(wfn.grid.directions);
    Real dX = wfn.grid.directions[0].getIncrementVector()[0];
    Real dY = wfn.grid.directions[1].getIncrementVector()[1];
    Real dZ = wfn.grid.directions[2].getIncrementVector()[2];
    Real dKX = 2 * M_PI * dZ / (dX * Real(nX));
    Real dKY = 2 * M_PI * dZ / (dY * Real(nY));
    Real E = wfn.getEnergy() * dZ * dZ;

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
        planeFourier /= nX*nY;

        planeFourier *= nX*nY;
        std::vector<Real> fourierReal(nXF*nYF);
        for(int i =0;i<nXF;++i)
        for(int j =0;j<nYF;++j)
            fourierReal[i*nYF+j] = abs(planeFourier(i,j));
        // Plotting the coefficients2 of the hartley transform
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
    else if ( mode == isoSurface ) {
        
        int nKX = nX, nKY = nY;
        int nK  = nKX * nKY;
        Real kXMax = nKX/2.0 * dKX;
        Real kYMax = nKY/2.0 * dKY;
        
        // Take only low-freq. Fourier components such that no basis function
        // decays faster than f(z) \propto 10^{-minExpGoal * z}
        std::cout << "Original wave vector grid: kx,ky = " 
                  << nX << "," << nY << "\n";
        Uint zSpan = zSurfEndIndex - zStartIndex;
        
        // Translate decayCutoff to minimum exponent per z-index
        Real minExpGoal = this->decayCutoff * dZ * log(10.0);   
        Real scaling = sqrt( (minExpGoal * minExpGoal + 2*E)
                             / (kXMax * kXMax + kYMax * kYMax) );
        if(scaling < 1.0){
            nKX = int( scaling * nX); kXMax = nKX/2.0 * dKX;
            nKY = int( scaling * nY); kYMax = nKY/2.0 * dKY;
            nK = nKX * nKY;
            std::cout << "Reduced wave vector grid: kx,ky = " << nKX << "," << nKY << "\n";
        }
        Real minExp = - sqrt( -2.0 * E + kXMax * kXMax + kYMax * kYMax);
        std::cout << "Strongest considered decay: x " << exp(minExp/dZ) << " per a.u.\n";
        minExp *= Real(zSurfEndIndex) - Real(zStartIndex);
        std::cout << "Smallest exponential term in matrix: " << exp(minExp) << "\n";
        
        
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
        //              $                   WORK, LWORK, IWORK, INFO )
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
//        String hartleyFile = io::getFileName(wfn.getFileName()) + ".hartley";
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

//        // Plotting the coefficients2 of the hartley transform
//        s = formats::gnuplot::writeMatrix<Real>(hartley, nX, nY);
//        hartleyFile = io::getFileName(wfn.getFileName()) + ".hartley2";
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

//        std::vector<Real> fourierReal(fourier.size());
//        for(int i =0;i<fourierReal.size();++i)
//            fourierReal[i] = abs(fourier[i]);
//        // Plotting the coefficients2 of the hartley transform
//        s = formats::gnuplot::writeMatrix<Real>(fourierReal, nXF, nYF);
//        hartleyFile = io::getFileName(wfn.getFileName()) + ".fourier";
//        io::writeStream(hartleyFile, s);

        
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

            // The c2r transform destroys its input (also the hartley?? ->check)
            tempFourier = planeFourier;
            plan_backward = fftw_plan_dft_c2r_2d(
                    nX, nY, 
                    (fftw_complex*) tempFourier.data(), 
                    tempDirect.data(), 
                    FFTW_ESTIMATE);

            // Perform Fourier transform backwards
            fftw_execute(plan_backward);

            //std::vector<Real> tmp; tmp.assign(tempDirect.begin(), tempDirect.end());
            //tempDirect /= nX*nY / (nKX*nKY);

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


    std::cout << "Time to extrapolate : " << (clock() -t)/1000.0 << " ms\n";

}

void WfnExtrapolation::determineRange(){
    using namespace blitz;

    WfnCube tmp(wfn);

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
        this->zSurfEndIndex = indices[2];
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
        this->zIndices = std::vector<Uint> (wfn.nX() * wfn.nY(), zStartIndex);
    }

    else if ( mode == isoSurface ) {
        if (approachFrom < 0) const_cast<Cube *>(hartree)->grid.zIsoSurfaceOnGrid(isoValue, this->zIndices);
        else                  const_cast<Cube *>(hartree)->grid.zIsoSurfaceOnGrid(isoValue, this->zIndices, approachFrom);
        wfn.grid.zSurface(this->zIndices, this->surface);
        
        types::String s = formats::gnuplot::writeMatrix<Real>(this->surface, wfn.nX(), wfn.nY());
        String surfaFile = io::getFileName(wfn.getFileName()) + ".exsurf";
        io::writeStream(surfaFile, s);
        s = formats::gnuplot::writeMatrix<Uint>(this->zIndices, wfn.nX(), wfn.nY());
        surfaFile = io::getFileName(wfn.getFileName()) + ".zindices";
        io::writeStream(surfaFile, s);
        
        Array<types::Uint,2> zIndicesBlitz(
                &zIndices[0],
                shape(wfn.nX(), wfn.nY()),
                neverDeleteData);
        
        zStartIndex = min(zIndicesBlitz);
        Real dZ = wfn.grid.directions[2].getIncrementVector()[2];
        zSurfEndIndex = max(zIndicesBlitz);
        zEndIndex = zSurfEndIndex + Uint(zWidth / dZ);

        std::cout << "Isosurface ranges over " << (zSurfEndIndex-zStartIndex)*dZ
                  << " from z = " << zStartIndex * dZ 
                  << " to " << max(zIndicesBlitz) * dZ << " [a.u.].\n";
        std::cout << "Extrapolation will be extended to "
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
   
    Cube* hptr = const_cast<Cube *>(hartree);

    std::vector<Real> hartreeVector;
    hptr->grid.zSurface(zIndices, hartreeVector);
    Array<Real, 2> hartreeSurface (
            &hartreeVector[0],
            shape(hartree->nX(), hartree->nY()),
            neverDeleteData);
    Real surfacePotential = mean(hartreeSurface);
    Array<Real, 2> hartreeeGrid (
            &(hptr->grid[0]),
            shape(hartree->nX(), hartree->nY(), hartree-nZ()),
            neverDeleteData);
    Real vacuumPotential = min(hartreeGrid);

    std::cout << "Hartree potential on extrapolation surface is " 
              << surfacePotential << " Ha\n";
    std::cout << "This is " << vacuumPotential - surfacePotential << 
              << " Ha below the maximum Hartree potential.\n";

 
    wfn.setEnergy( wfn.getEnergy()- vacuumPotential);
    std::cout << "Energy level " << wfn.getEnergy() << " Ha\n";


    // Print hartree potential for gnuplot
    types::String s = formats::gnuplot::writeMatrix(hartreeVector, hartree->nX(), hartree->nY());
    String surfaceFile = io::getFileName(hartree->getFileName()) + ".exsurf";
    io::writeStream(surfaceFile, s);

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
