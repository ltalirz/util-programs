#include <iostream>
#include <string>

#include "atomistic.hpp"

using namespace std;

int main(){
    atomistic::Spectrum cp = atomistic::Spectrum();
    cp.readFromCp2k("stm.out");
    cp.print();

    return 0;
}

