#include "boost/progress.hpp"
#include<stdio.h>



int main() {
    boost::progress_display d(100);
    for(int i = 0; i<100; ++i){
        sleep(1);
        std::cout << "i\n";
        ++d;
    }
    return 0;
}
