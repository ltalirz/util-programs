#include <iostream>
#include <string>

#include "types.hpp"

using namespace std;

int main(){
    throw types::runtimeError() << types::errinfo_runtime("Cannot read cube file - no file name given.");
    throw types::parseError() << types::errinfo_parse("title or description");

    return 0;
}

