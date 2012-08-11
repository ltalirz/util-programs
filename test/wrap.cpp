// Testing the dge_pseudo wrapper that was written
// for extrapolate2 in order to perform an efficient
// pseudoinverse.

#include<iostream>
#include<vector>
#include<math.h>

#include "wrappers/lapack.hpp"

using namespace wrappers::lapack;

void printMatrix(std::vector<double>& m, int M, int N){
    for(int i=0; i<M; ++i){
        for(int j=0; j<N; ++j){
            std::cout << m[i*N + j] << " ";
        }
        std::cout << std::endl;
    }
}


void transposeMatrix(std::vector<double>& m, int M, int N){
    std::vector<double> tmp = m;
    for(int i=0; i<M; ++i){
        for(int j=0; j<N; ++j){
            m[i*N+j] = tmp[j*M+i];
        }
    }
}


int main(){
    int M=4, N=5;
    std::vector<double> m(M*N, 0.0);
    m[0] = 1;
    m[4] = 2;
    m[N+2] = 3;
    m[3*N+1] = 4;

    printMatrix(m, M, N);

    transposeMatrix(m, M, N);
    printMatrix(m, N, M);
    dge_pseudo(m, M, N, 1e-5);
    printMatrix(m, M, N);

    return 0;
}


