// Testing the dge_pseudo wrapper that was written
// for extrapolate2 in order to perform an efficient
// pseudoinverse.

// See the Wikipadia article on single value decomposition
// for the values used in this example.

#include<iostream>
#include<vector>
#include<math.h>

#include "wrappers/lapack.hpp"

using namespace wrappers::lapack;

void printMatrixRM(std::vector<double>& m, int M, int N){
    for(int i=0; i<M; ++i){
        for(int j=0; j<N; ++j){
            std::cout << m[i*N + j] << " ";
        }
        std::cout << std::endl;
    }
}

void printMatrixCM(std::vector<double>& m, int M, int N){
    for(int i=0; i<M; ++i){
        for(int j=0; j<N; ++j){
            std::cout << m[j*M + i] << " ";
        }
        std::cout << std::endl;
    }
}


void transposeMatrix(std::vector<double>& m, int M, int N){
    std::vector<double> tmp = m;
    for(int i=0; i<M; ++i){
        for(int j=0; j<N; ++j){
            m[j*M+i] = tmp[i*N+j];
        }
    }
}


int main(){

    // See example on wikipedia
    int M=4, N=5;
    std::vector<double> m(M*N, 0.0);
    m[0] = 1;
    m[4] = 2;
    m[N+2] = 3;
    m[3*N+1] = 4;

    std::cout << "A\n";
    printMatrixRM(m, M, N);

    bool moreColumnsThanRows = true;

    if(moreColumnsThanRows){
        // Taking the matrix as it is
        // (the transpose is needed for the fortran
        // columnn-major format)

        transposeMatrix(m, M, N);
        dge_pseudo(m, M, N, 1e-5);
        std::cout << "Pseudoinverse of A\n";
        printMatrixCM(m, N, M);
    }
    else{
        // Feeding it the transposed matrix.
        // (since we have stored the matrix in row-major format,
        // we don't need a transposition here)

        dge_pseudo(m, N, M, 1e-5);
        std::cout << "Pseudoinverse of A^T\n";
        printMatrixCM(m, N, M);
    }
    return 0;
}


