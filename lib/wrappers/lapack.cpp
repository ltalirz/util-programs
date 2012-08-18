#include "wrappers/lapack.hpp"
#include "types.hpp"

namespace wrappers {
namespace lapack {


int dge_pseudo(std::vector<double> &A, int M, int N, double RCOND){
    
    char JOBU;
    char JOBVT;
    int  LDA  = M;
    std::vector<double> S(0);
    std::vector<double> U(0);
    int  LDU  = M;
    std::vector<double> VT(0);
    int  LDVT  = N;
    std::vector<double> WORK(1);
    int  LWORK  = -1;
    int  INFO   = 0;
    
    if (M >= N){
        JOBU = 'O';
        JOBVT = 'A';
        VT = std::vector<double>(N*N);
        S = std::vector<double>(N);
    }
    else{
        JOBU = 'A';
        JOBVT = 'O';
        U  = std::vector<double>(M*M);
        S = std::vector<double>(M);
    }

    // Find out optimum WORK size
    dgesvd_(&JOBU, &JOBVT,
            &M, &N,
            &*A.begin(), &LDA,
            &*S.begin(),
            &*U.begin(), &LDU,
            &*VT.begin(), &LDVT,
            &*WORK.begin(), &LWORK,
            &INFO);
    
    LWORK = int(abs(WORK[0])); 
    WORK =  std::vector<double>(LWORK);
    
    // Now perform A = U * S * VT
    dgesvd_(&JOBU, &JOBVT,
            &M, &N,
            &*A.begin(), &LDA,
            &*S.begin(),
            &*U.begin(), &LDU,
            &*VT.begin(), &LDVT,
            &*WORK.begin(), &LWORK,
            &INFO);

    if (INFO != 0) return INFO;

    // A+ = V * S+ * UT
    // First calculate S+, setting small singular values to zero
    std::vector<double>::iterator it = S.begin(), end = S.end();
    double SMAX = *it;
    int    NS   = 0;
    while(it != end){
        if(*it / SMAX < RCOND) break;
        else *it = 1.0 / *it;
        ++it;
        ++NS;
    }
    while(it != end){
        *it = 0.0;
        ++it;
    }
    
    // Now multiply S+ to the smaller matrix and then
    // the result to the larger matrix.

    std::vector<double> C(N*M);
    if (M >= N){
        // Need V * S+ (general * diagonal), which is not a BLAS.
        // Since VT is column-major, we can omit transposition
        // and just think the C++ way.
        for(int i = 0; i < NS; ++i)
            for(int j = 0; j < N; ++j)
                VT[j * N + i] *= S[i];

        // The result is nominally a NxM matrix, but actually only
        // NS columns are nonzero (and N columns stored)

        // A+ = V * S+ * UT
        // A stores the first N columns of U, 
        // i.e. AT has the first N rows of UT
        char TRANSA = 'T';
        char TRANSB = 'T';
        int K  = NS;
        double ALPHA = 1.0;
        int LDA = N;
        int LDB = M;
        double BETA = 0.0;
        int LDC = N;

        dgemm(&TRANSA, &TRANSB,
                &N,&M,&K,&ALPHA,
                &*VT.begin(),&LDA,
                &*A.begin(),&LDB,
                &BETA,&*C.begin(),&LDC);
    }
    else{
        // Need S+ * UT (diagonal * general), which is not a BLAS.
        // Since U is column-major, we can omit transposition
        // and just think the C++ way.
        for(int i = 0; i < NS; ++i)
            for(int j = 0; j < M; ++j)
                U[i*M + j] *= S[i];

        // The result is nominally a NxM matrix, but actually only
        // NS rows are nonzero (and N rows stored)

        // A+ = V * S+ * UT
        // A stores the first M columns of VT
        // i.e. AT has the first M rows of V
        char TRANSA = 'T';
        char TRANSB = 'T';
        int K  = NS;
        double ALPHA = 1.0;
        int LDA = M;  // Note: In column-major the leading dimension is the
        int LDB = M;  //       number of rows. Need ld *before* transposition.
        double BETA = 0.0;
        int LDC = N;

        dgemm(&TRANSA, &TRANSB,
                &N,&M,&K,&ALPHA,
                &*A.begin(),&LDA,
                &*U.begin(),&LDB,
                &BETA,&*C.begin(),&LDC);
        
    }
    
    // Overwrite A with its pseudoinverse
    A = C;
    
    return 0;


}



}
}
