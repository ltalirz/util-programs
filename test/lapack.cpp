#include<iostream>
#include<complex>
#include<math.h>

#define __MKL__

#ifdef __MKL__
#define MKL_Complex16 std::complex<double>
#include "mkl.h"

#else
extern "C" void zgetrf_( int*, int* , complex<double>* , int*, int* , int* );
extern "C" void zgetri_( int*, complex<double>* , int*, int* , complex<double>*, int* , int* );

#endif

using namespace std;
typedef std::complex<double> cmplx;

//........................................................................................
void zgeTranspose( cmplx *Transposed, cmplx *M ,int n)
{

    int i,j;
    for(i=0;i<n;i++)

        for(j=0;j<n;j++) Transposed[i+n*j] = M[i*n+j];
}

//.........................................................................................
void MatrixComplexInverse(cmplx *invA, cmplx *A, int n)
{

    int LWORK=10*n;

    int *permutations;

    cmplx *WORK, *tempA;

    tempA = new cmplx[n*n];

    permutations =  new int[2*n];
    WORK = new cmplx[n*n];


    int INFO;

    zgeTranspose(tempA,A,n);


    zgetrf_( &n, &n, tempA , &n, permutations , &INFO );

    if (INFO != 0) {
        cout<<"ComplexMatrixInverse: Error at zgetrf  \n"; exit(0);
    }



    zgetri_( &n, tempA , &n, permutations , WORK, &LWORK, &INFO );

    if (INFO != 0) {
        cout<<"ComplexMatrixInverse: Error at zgetri  \n"; exit(0);
    }

    zgeTranspose(invA,tempA,n);

    delete [] WORK;

    delete [] tempA;
    delete [] permutations;

}

/////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    int i,j;
    const int N = 3;

    cmplx I(0.,1.);

    cmplx A[] = { 1. + I , 2. ,  3 , 4. , 5.+I , 6. , 7., 8., 9. + I};

    cmplx invA[N*N];

    MatrixComplexInverse(invA,A,N);


    for(i=0;i<N;i++){
        for(j=0;j<N;j++) cout << invA[i*N + j]<<"\t";

        std::cout<<"\n";
    }

    std::cout<<"---------------------------\n";

    return 0;

}
