/*
 * gnuplot.h
 */
#ifndef FORMATS_GNUPLOT_H
#define FORMATS_GNUPLOT_H

#include "types.hpp"
#include <boost/format.hpp>

namespace formats {

namespace gnuplot {

template<typename T> types::String writeMatrix(std::vector<T> data, types::Uint nX, types::Uint nY){
    typename std::vector<T>::const_iterator dataIt = data.begin();
    types::String result = "";
    for(types::Uint x = 0; x<nX; ++x){
        for(types::Uint y=0; y<nY; ++y){
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

#endif
