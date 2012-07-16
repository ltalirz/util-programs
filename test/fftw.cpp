#include <complex>
#include <fftw3.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <vector>

// Handcoded 1d DFT on a 2d function
void r2c1dz_code(){
    size_t N = 100;

    // Here we generate a 2d function that decays in z like exp(-kz z)
    // where kz = sqrt(2E - kx^2)
    // We choose E = 5/100 and kx = 1/10, yielding kz = 2/10 
    typedef std::complex<double> cmplx; 

    
    // x (j) is the fast index, z (i) the slow one
    // WE CHOOSE KX = 1/n
    // => want kz=2/N => E = 5/2 /N/N
    double kx = 1.0/N;
    double E  = 1.5/N/N;
    double kz = sqrt(2 * E + kx*kx); // = 2/N
    std::cout << kz;
    
    std::vector< double > input;
    for(int z = 0; z < N; ++z) {
        for(int x = 0; x < N; ++x){
            input.push_back( sin(2*3.141*x * kx)*exp(-kz*z)  );
            fprintf( stdout, "data[%d, %d] = { %2.2f }\n",
                 z, x, sin(2*3.141*x/(double) N)*exp(-kz*z) );
        }
    }
   

    std::vector< cmplx > dft;
    cmplx          temp  = 0;

    std::vector< int > st;
    for(int i = 0; i < N; ++i) {
         // stutzpunlte alle glecih bei z=0
         st.push_back(i);
    }

    // Perform DFT
       double z;
    for(int kxj = 0; kxj < N; ++kxj){
       temp = 0;
       if (kxj >= N/2) kx = (double)(N-kxj)/N;
       else            kx = (double)kxj/N;
       double mykz = sqrt( 2*E + kx * kx );
       std::cout << mykz;
        for(int x = 0; x < N; ++x){
            z = double(st[x]);
            temp += input[x + N*z] * exp( 2 * M_PI * cmplx(0,1) * double(x) * double(kxj)/double(N) + mykz * z);
            //std::cout <<  input[x + N*z] * exp( 2 * M_PI * cmplx(0,1) * double(x) * double(kx)/double(N) + mykz * z);
            //std::cout <<  input[x] * exp( 2 * M_PI * cmplx(0,1) * double(x) * double(kx)/double(N) );
        }
        dft.push_back(temp);
    }

    /* print dft result */
    for(int x = 0 ; x < N ; x++ ) {
            fprintf( stdout, "fft_result[%d] = { %2.2f }\n",
                 x, std::abs(dft[x]) );
    }

}


// Handcoded 1d DFT
void r2c1d_code(){
    size_t N = 10;

    typedef std::complex<double> cmplx; 
//
    std::vector< double > input;
    for(int i = 0; i < N; ++i) {
            input.push_back( sin(2*3.141*i/(double) N) );
            fprintf( stdout, "data[%d] = { %2.2f }\n",
                 i, sin(2*3.141*i/(double) N) );
     }
    
    std::vector< cmplx > dft;
    cmplx          temp  = 0;

    // Perform DFT
    for(int i = 0; i < N; ++i){
        temp = 0;
        for(int j = 0; j < N; ++j){
            temp += input[j] * exp(-2 * M_PI * cmplx(0,1) * double(j) * double(i)/double(N));
        }
        dft.push_back(temp);
    }

    /* print dft result */
    for(int i = 0 ; i < N ; i++ ) {
            fprintf( stdout, "fft_result[%d] = { %2.2f }\n",
                 i, std::abs(dft[i]) );
    }
    
    double *in = &input[0];
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2 + 1));

    fftw_plan plan_forward, plan_backward;

    plan_forward = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);
    fftw_execute(plan_forward); /* repeat as needed */


    fftw_execute( plan_forward );

    /* print fft result */
    std::complex<double> *pt = reinterpret_cast< std::complex<double>* >(out);
    for(int i = 0 ; i < N ; i++ ) {
            fprintf( stdout, "fft_result[%d] = { %2.2f }\n",
                 i,  std::abs(*pt) );
            ++pt;
    }

}



void r2c2d() {
    size_t N = 10;

    std::vector< double > input;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j){
            input.push_back( sin(2*3.141*i/(double) N)*sin(2*3.141*j/(double) N)  );
            fprintf( stdout, "data[%d, %d] = { %2.2f }\n",
                 i, j, sin(2*3.141*i/(double) N)*sin(2*3.141*j/(double) N) );
        }
    }
    double *in = &input[0];
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * (N/2 + 1));

    fftw_plan plan_forward, plan_backward;

    plan_forward = fftw_plan_dft_r2c_2d(N,N, in, out, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_c2r_2d(N,N, out, in, FFTW_ESTIMATE);
    fftw_execute(plan_forward); /* repeat as needed */


    fftw_execute( plan_forward );

    /* print fft result */
    std::complex<double> *pt = reinterpret_cast< std::complex<double>* >(out);
    for(int i = 0 ; i < N ; i++ ) {
        for(int j = 0; j < N/2 +1; ++j){
            fprintf( stdout, "fft_result[%d,%d] = { %2.2f }\n",
                 i, j, std::abs(*pt) );
            ++pt;
        }
    }

//    fftw_execute( plan_backward );
//
//    /* print ifft result */
//    for(int i = 0 ; i < N ; i++ ) {
//        fprintf( stdout, "ifft_result[%d] = { %2.2f, %2.2f }\n",
//                 i, in[i][0] / N, in[i][1] );
//    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
//    fftw_free(in);
    fftw_free(out);
}



// Perform fft on top of STL C++ containers
void usingStdTypes() {
    size_t N = 10;

    std::vector< std::complex<double> > input;
    for(int i = 0; i < N; ++i) {
        input.push_back( std::complex<double>( sin(2*3.141*i/(double) N), 0.0 ) );
        fprintf( stdout, "data[%d] = { %2.2f, %2.2f }\n",
                 i, input[i].real(), input[i].imag() );
    }
    fftw_complex *in = reinterpret_cast<fftw_complex*>(&input[0]);
    fftw_complex *out;

    fftw_plan plan_forward, plan_backward;

    out =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    plan_forward = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward); /* repeat as needed */


    fftw_execute( plan_forward );

    /* print fft result */
    for(int i = 0 ; i < N ; i++ ) {
        fprintf( stdout, "fft_result[%d] = { %2.2f, %2.2f }\n",
                 i, out[i][0], out[i][1] );
    }

    fftw_execute( plan_backward );

    /* print ifft result */
    for(int i = 0 ; i < N ; i++ ) {
        fprintf( stdout, "ifft_result[%d] = { %2.2f, %2.2f }\n",
                 i, in[i][0] / N, in[i][1] );
    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
//    fftw_free(in);
    fftw_free(out);
}

// Perform fft on top of fftw-specialized containers
void usingFftwTypes() {
    fftw_complex *in, *out;
    fftw_plan plan_forward, plan_backward;
    size_t N = 10;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    for(int i = 0; i < N; ++i) {
        in[i][0] = sin(2*3.141*i/(double) N);
        in[i][1] = 0.0;
        fprintf( stdout, "data[%d] = { %2.2f, %2.2f }\n",
                 i, in[i][0], in[i][1] );
    }

    out =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    plan_forward = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_forward); /* repeat as needed */


    fftw_execute( plan_forward );

    /* print fft result */
    for(int i = 0 ; i < N ; i++ ) {
        fprintf( stdout, "fft_result[%d] = { %2.2f, %2.2f }\n",
                 i, out[i][0], out[i][1] );
    }

    fftw_execute( plan_backward );

    /* print ifft result */
    for(int i = 0 ; i < N ; i++ ) {
        fprintf( stdout, "ifft_result[%d] = { %2.2f, %2.2f }\n",
                 i, in[i][0] / N, in[i][1] );
    }

    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_free(in);
    fftw_free(out);
}

int main() {
//    usingStdTypes();
//    usingFftwTypes();
//  r2c2d();
//   r2c1d_code();
    r2c1dz_code();
    return 0;
}
