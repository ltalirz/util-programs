#include "io.h"
#include "atomistic.h"
#include "types.h"

#include <iostream>
#include <vector>
#include <ctime>

using namespace std;
using namespace types;

const string filenameOut = "out.cube";
const string filename = "2MOL-WFN_00122_1-1_23.cube";
const int NFLOATS = 47000000; // The number of floats in the big metal+mol cubefile

// Test String concatenation (try to get below 30s on this)
// Reserving string size does not help
// Operating on std::vector<char> instead of std::string doesn'OBt either
void addData() {

    using boost::spirit::karma::right_align;
    using boost::spirit::karma::repeat;
    using boost::spirit::karma::int_;
    using boost::spirit::karma::columns;
    using boost::spirit::karma::eol;
    using boost::spirit::karma::generate;

    Binary stream; stream.reserve(47000000 * 20);
    std::back_insert_iterator<Binary> sink(stream);
    std::vector<double> data; data.reserve(NFLOATS);
    for(int i=0; i<NFLOATS;++i){
            data.push_back(i*3);
    }
    
    std::vector<Real>::const_iterator dataIt = data.begin();
    // Fastest direction is z, stored in directions[2]
    Uint i = 1, nZ = 20;
    time_t t = clock();
    while(dataIt != data.end()) {
        generate(sink, right_align(13)[types::sci5] << right_align(13)[types::sci5],  *dataIt, *(++dataIt));
//        if(i % 6 == 0) generate(sink, eol);
//        else if(i % nZ == 0) {
//            generate(sink, eol);
//            i=0;
//        }
        ++dataIt;
        ++i;
    }
    std::cout << "Time to add : " << (clock() -t)/1000.0 << " ms\n";

}






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
//  readWriteCubeFile();
    addData();
    return 0;
}

