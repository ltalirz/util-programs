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
    
    if (M > N){
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

    std::cout << "Performed SVD\n";
    
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
    
    std::cout << NS << "Inverted Sigma\n";

    // Then multiply S+ to the smaller matrix
    // Since diagonal * general matrix is not a BLAS,
    // we implement it ourselves.
    std::vector<double> C(N*M);

    if (M > N){
        // Need V * S+
        // Since VT is column-major, we are automatically
        // transposing it.
        for(int i = 0; i < NS; ++i)
            for(int j = 0; j < N; ++j)
                VT[j * N + i] *= S[i];

        // The result is nominally a NxM matrix, but actually only
        // NS columns are nonzero (and N columns stored)

        // A+ = V * S+ * UT
        char TRANSA = 'N';
        char TRANSB = 'T';
        int K  = NS;
        double ALPHA = 1.0;
        int LDA = N;
        int LDB = M;
        double BETA = 0.0;
        int LDC = N;

        dgemm(&TRANSA, &TRANSB,
                &M,&N,&K,&ALPHA,
                &*VT.begin(),&LDA,
                &*U.begin(),&LDB,
                &BETA,&*C.begin(),&LDC);
    }
    else{
        // Need S+ * UT
        // Since UT is column-major, we are automatically
        // transposing it by thinking the C++ way.
        for(int i = 0; i < NS; ++i)
            for(int j = 0; j < N; ++j)
                U[j * N + i] *= S[i];

        // The result is nominally a NxM matrix, but actually only
        // NS rows are nonzero (and N rows stored)

        // A+ = V * S+ * UT
        char TRANSA = 'T';
        char TRANSB = 'N';
        int K  = NS;
        double ALPHA = 1.0;
        int LDA = N;
        int LDB = N;
        double BETA = 0.0;
        int LDC = N;

        dgemm(&TRANSA, &TRANSB,
                &M,&N,&K,&ALPHA,
                &*U.begin(),&LDA,
                &*VT.begin(),&LDB,
                &BETA,&*C.begin(),&LDC);
    }

    // Overwrite A with its pseudoinverse
    A = C;
    
    return 0;


}



}
}
