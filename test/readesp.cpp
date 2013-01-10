#include <iostream>
#include <string>

#include "formats.hpp"

using namespace std;

int main(){
    formats::espresso::Spectrum s = formats::espresso::Spectrum();
    s.readFromOutput("esp.out");
    s.print();

    return 0;
}

