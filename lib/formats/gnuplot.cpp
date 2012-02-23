#include "formats/gnuplot.h"
#include "types.hpp"
#include <boost/format.hpp>


namespace formats {
    namespace gnuplot {

        using namespace types;

        String writeMatrix(std::vector<Real> data, Uint nX, Uint nY){
            std::vector<Real>::const_iterator dataIt = data.begin();
            String result = "";
            for(Uint x = 0; x<nX; ++x){
                for(Uint y=0; y<nY; ++y){
                    result += str(boost::format("%12.6e") % *dataIt);
                    result += " ";
                    ++dataIt;
                }
                result += "\n";
            }
            return result;
        }


    }
}
