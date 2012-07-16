#include <complex>
#include <fftw3.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <blitz/array.h>
#include "lapack.hpp"


typedef std::complex<double> cmplx; 


void printMatrix(const std::vector<cmplx>& matrix){
    size_t n = sqrt(matrix.size());
    std::vector<cmplx>::const_iterator it = matrix.begin();
    for(size_t i = 0; i < n; ++i){
        std::cout << " | ";
        for(size_t j = 0; j < n; ++j){
            fprintf( stdout, " %+2.2f %+2.2f i ", real(*it), imag(*it));
            ++it;
        }
        std::cout << " | \n";
    }
    std::cout << std::endl;
}
void printVector(const std::vector<cmplx>& line){
    size_t n = line.size();
    std::vector<cmplx>::const_iterator it = line.begin();
        std::cout << " | ";
        for(size_t j = 0; j < n; ++j){
            fprintf( stdout, " %+2.2f %+2.2f i ", real(*it), imag(*it));
            ++it;
        }
        std::cout << " | \n";
    std::cout << std::endl;
}

// Handcoded 1d DFT on a 2d function
void r2c1dz_code_v2(){
    size_t N = 10000;
    // v1: N = 10 000 needed 7.5 GB RES MEM and 15 minutes (infinijazz node)
    
    // Here we generate a 2d function input(x,z) = input[x+z*N]
    // that decays in z like exp(-kz z) where kz = sqrt(2E + kx^2)   
    // We CHOOSE kx = 1/N and want kz=2/N => E = 3/2 /N/N
    // Goal: Obtain the Fourier coefficients of this function in x,
    //  no matter at which z we give the function.
    double kx = 1.0/N;
    double E  = 1.5/N/N;
    double kz = sqrt(2 * E + kx*kx); // = 2/N

    std::vector< cmplx > input;
    for(int z = 0; z < N; ++z) {
        for(int x = 0; x < N; ++x){
            input.push_back( sin(2.0 * M_PI *x * kx)*exp(-kz*z)  );
        }
    }
    //printMatrix(input); 


    std::vector< int > st (N);   // z-indices, where we know f
    std::vector<cmplx> f  (N);   // the function value f at these points

    for(int i = 0; i < N; ++i) {
         st[i] = i;              // lets say we know f on a "stair" 
         f[i]  = input[ i + N* st[i]];
    }
    
    // Build matrix for modified Fourier transform
    // now need transposed matrix
    std::vector<cmplx> matrix (N*N);
    for(int i =0; i < N; ++i)
        for(int j=0; j < N; ++j){
            double k = (2*j < N)? double(j)/N : (double(N) - j)/N;
            matrix[i + j*N] = 
                exp( 2.0 * M_PI * cmplx(0,1) * double(i) * double(j)/double(N) 
                     - sqrt( 2*E + k*k) * st[i]);

        }
    //printMatrix(matrix);

    // Invert this fucker
    int one=1;
    int info;
    int n = N;
    std::vector< int > ipiv (N);   // z-indices, where we know f
    std::cout << "Before inverse" << std::endl;
    zgesv_(&n, &one, &*matrix.begin(), &n, &*ipiv.begin(), &*f.begin(), &n, &info);

    //printVector(coeff); 
    std::cout << " coeff[1] " << f[1] << " coeff[N-1] " << f[N-1];

}    


// Handcoded 1d DFT on a 2d function
void r2c1dz_code(){
    size_t N = 10000;
    // v1: N = 10 000 needed 7.5 GB RES MEM and 15 minutes (infinijazz node)
    
    // Here we generate a 2d function input(x,z) = input[x+z*N]
    // that decays in z like exp(-kz z) where kz = sqrt(2E + kx^2)   
    // We CHOOSE kx = 1/N and want kz=2/N => E = 3/2 /N/N
    // Goal: Obtain the Fourier coefficients of this function in x,
    //  no matter at which z we give the function.
    double kx = 1.0/N;
    double E  = 1.5/N/N;
    double kz = sqrt(2 * E + kx*kx); // = 2/N

    std::vector< cmplx > input;
    for(int z = 0; z < N; ++z) {
        for(int x = 0; x < N; ++x){
            input.push_back( sin(2.0 * M_PI *x * kx)*exp(-kz*z)  );
        }
    }
    //printMatrix(input); 


    std::vector< int > st (N);   // z-indices, where we know f
    std::vector<cmplx> f  (N);   // the function value f at these points

    for(int i = 0; i < N; ++i) {
         st[i] = i;              // lets say we know f on a "stair" 
         f[i]  = input[ i + N* st[i]];
    }
    
    // Build matrix for modified Fourier transform
    std::vector<cmplx> matrix (N*N);
    for(int i =0; i < N; ++i)
        for(int j=0; j < N; ++j){
            double k = (2*j < N)? double(j)/N : (double(N) - j)/N;
            matrix[i*N+j] = 
                exp( 2.0 * M_PI * cmplx(0,1) * double(i) * double(j)/double(N) 
                     - sqrt( 2*E + k*k) * st[i]);

        }
    //printMatrix(matrix);

    // Invert this fucker
    std::vector<cmplx> inverse(N*N);
    std::cout << "Before inverse" << std::endl;
    MatrixComplexInverse(&inverse[0],&matrix[0],  N);     
    //printMatrix(inverse);
    std::cout << "afterinverse";

    // Obtain Fourier coefficients
    std::vector<cmplx> coeff(N, 0.0);

    for (int i = 0; i<N; ++i){
        for(int j = 0; j < N; ++j){
            coeff[i] += inverse[i*N + j] * f[j];
        }
    }
    //printVector(coeff); 
    std::cout << " coeff[1] " << coeff[1] << " coeff[N-1] " << coeff[N-1];

}    
    
// Handcoded 2d DFT on a 3d function
void r2c2dz_code(){
    size_t dimX = 10;
    size_t dimY = 10;
    size_t N = dimX * dimY;
     
    // Here we generate a 3d function input(x,y,z) = input[x+y*N+z*N*N]
    // that decays in z like exp(-kz z) where kz = sqrt(2E + kx^2 + ky^2)   
    // We CHOOSE kx = 1/N = ky and want kz=2/N => E = 1/N/N
    // Goal: Obtain the Fourier coefficients of this function in x,
    //  no matter at which z we give the function.
    double kx = 1.0/dimX;
    double ky = 1.0/dimY;
    double E  = 1.0/N;
    double kz = sqrt(2 * E + kx*kx + ky*ky); // = 2/N

    std::vector< cmplx > input;
    for(int z = 0; z < N; ++z) {
        for(int y = 0; y < dimY; ++y) {
            for(int x = 0; x < dimX; ++x){
                input.push_back( sin(2.0 * M_PI *(x * kx+ y*ky)) * exp(-kz*z) );
            }
        }
    }
    //printMatrix(input); 


    std::vector< int > st (N);   // z-indices, where we know f
    std::vector<cmplx> f  (N);   // the function value f at these points

    for(int y = 0; y < N; ++y) {
        for(int i = 0; i < N; ++i) {
         st[i] = i;              // lets say we know f on a "stair" 
         f[i]  = input[ i + N* st[i]];
    }
    
    // Build matrix for modified Fourier transform
    std::vector<cmplx> matrix (N*N);
    for(int i =0; i < N; ++i)
        for(int j=0; j < N; ++j){
            double k = (2*j < N)? double(j)/N : (double(N) - j)/N;
            matrix[i*N+j] = 
                exp( 2.0 * M_PI * cmplx(0,1) * double(i) * double(j)/double(N) 
                     - sqrt( 2*E + k*k) * st[i]);

        }
    printMatrix(matrix);

//    // Invert this fucker
//    std::vector<cmplx> inverse(N*N);
//    MatrixComplexInverse(&inverse[0],&matrix[0],  N);     
//    //printMatrix(inverse);
//
//    // Obtain Fourier coefficients
//    std::vector<cmplx> coeff(N, 0.0);
//
//    for (int i = 0; i<N; ++i){
//        for(int j = 0; j < N; ++j){
//            coeff[i] += inverse[i*N + j] * f[j];
//        }
//    }
//    printVector(coeff); 
//
}    

int main() {
//  r2c1dz_code();
    r2c1dz_code_v2();
    return 0;
}
