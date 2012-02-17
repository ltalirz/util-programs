/*
 * cube.h
 */
#ifndef FORMATS_STM_H
#define FORMATS_STM_H

#include "types.hpp"
#include "formats/cube.h"

namespace formats {

/**
 * STM image
 * List of WfnCubes and bias
 */
class StmImage {
    private:
        std::vector<WfnCube> cubes;
        types::Real bias;

    public:
        StmImage(types::Real b) : bias(b) {};
        StmImage() {};
        bool findCubes(std::vector<WfnCube> &list);
        types::Real getBias() { return bias; }
};

}

#endif
