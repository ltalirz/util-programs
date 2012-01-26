#include "io.h"
#include "atomistic.h"
#include "types.h"

#include <iostream>
#include <vector>

using namespace std;

const string filenameOut = "out.cube";
const string filename = "2MOL-WFN_00122_1-1_23.cube";
const int NFLOATS = 3000000;

// Test the read/write CubeFile function from atomistic
// as well as resize 
void readWriteCubeFile() {
   atomistic::Cube cub = atomistic::Cube();
   cub.readCubeFile(filename);
    
   std::vector<types::Uint> increments;
   increments.push_back(10);
   increments.push_back(10);
   increments.push_back(100);

//cub.print();
   std::cout<< "Before first write \n";
   cub.writeCubeFile(filenameOut);
   std::cout<< "After  first write, before resize \n";
   
//   cub.grid.resize(increments); 
//   std::cout<< "After resize, before second write \n";
//   cub.writeCubeFile(filenameOut);
}


// Test the basic read/write function from io
void readWrite() {
    std::vector<char> data;
    io::readBinary(filename, data);
    io::writeBinary(filenameOut, data);
}


int main() {

//  readWrite();
    readWriteCubeFile();

    return 0;
}

