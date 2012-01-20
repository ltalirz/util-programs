#include <blitz/array.h>
#include <iostream>
#include <vector>

// Resizing array - works.
// Note: If array deletes data when done you get an error because of 
// double free (from vector and Array)
void resizeArray(unsigned int size) {
    using namespace blitz;
    
    std::vector<double> vec(size,3.4);
    
    Array<double,1> vecArray(&vec[0], shape(size), neverDeleteData);
   
    vecArray(1) = 100; 
    vecArray(0) = 70; 
    cout << vecArray << endl;
    vecArray.resizeAndPreserve(shape(size+1));
    cout << vecArray << endl;
}

// Partially mapping a vector works as well.
void doubleReuseVector(unsigned int size) {
    using namespace blitz;
    
    std::vector<double> vec(size,3.4);
    
    Array<double,1> vecArray(&vec[0], shape(size), neverDeleteData);
    Array<double,1> vecArray2(&vec[size/2], shape(size/2), neverDeleteData);
    vecArray2(1) = 20;
    vec[7]=100;
    cout << vecArray;
}


// This works perfectly. The Array reuses the vector. 
void reuseVector(unsigned int size) {
    using namespace blitz;
    
    std::vector<double> vec(size,3.4);
    Array<double,1> vecArray(&vec[0], shape(size), neverDeleteData);
    vec[7]=100;

    cout << vecArray;
}

void reuseArray(){
   using namespace::blitz; 
   using namespace tensor;
    double data[] = {1,2,3,4,5,6,7,8,9,10,11,12};    
    Array<double,3> dataArray(data, shape(3,2,2), neverDeleteData);
    cout << dataArray;
    std::cout << dataArray(0,0,1);
    cout << dataArray(0,1,0);
    cout << dataArray(1,0,0);

//    Array<double,1> reduced(sum(dataArray(j,i),j));
//    cout << reduced;
}

void reduceTensor() {
    using namespace blitz;
    Array<float,2> A(3,3);
    Array<float,1> B(3);
    A = 0, 1, 2,
    3, 4, 5,
    6, 7, 8;

    cout << A;
    B= sum(A(tensor::i,tensor::j),tensor::j);
    cout << B;
    
    cout << sum(A) << endl          // 36
         << min(A) << endl          // 0
         << count(A >= 4) << endl;  // 5
}

int main() {
    unsigned int size = 20;
//  doubleReuseVector(size);
//  reuseArray();
    resizeArray(size);
    return 0;
}
