/*
 * lapack.hpp
 */
#ifndef WRAPPERS_LAPACK_HPP
#define WRAPPERS_LAPACK_HPP

#include "types.hpp"

namespace wrappers {
namespace lapack   {


#define __MKL__

#ifdef __MKL__
#define MKL_Complex16 std::complex<double>
#include "mkl.h"

#else
extern "C" void zgetrf_( int*, int* , complex<double>* , int*, int* , int* );
extern "C" void zgetri_( int*, complex<double>* , int*, int* , complex<double>*, int* , int* );
#endif

// Calculates the pseudoinverse of M x N matrix A, 
// throwing away singular values smaller than RCOND * max(sv)
int dge_pseudo(std::vector<double> &A, int M, int N, double RCOND);

}
}

#endif
