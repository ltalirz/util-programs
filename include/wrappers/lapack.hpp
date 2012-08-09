/*
 * lapack.hpp
 */
#ifndef WRAPPERS_LAPACK_HPP
#define WRAPPERS_LAPACK_HPP

#define __MKL__

#ifdef __MKL__

#define MKL_Complex16 std::complex<double>
#include "mkl.h"

#else

extern "C" void zgetrf_( int*, int* , complex<double>* , int*, int* , int* );
extern "C" void zgetri_( int*, complex<double>* , int*, int* , complex<double>*, int* , int* );

#endif




#endif
