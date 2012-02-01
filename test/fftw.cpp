#include <complex>
#include <fftw3.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <vector>

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
    r2c2d();
    return 0;
}
