/*
 * gnuplot.h
 */
#ifndef FORMATS_GNUPLOT_H
#define FORMATS_GNUPLOT_H

#include "types.hpp"

namespace formats {

namespace gnuplot {
    types::String writeMatrix(
            std::vector<types::Real> data,
            types::Uint nX,
            types::Uint nY);
}

}

#endif
