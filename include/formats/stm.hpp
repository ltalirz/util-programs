/*
 * cube.h
 */
#ifndef FORMATS_STM_H
#define FORMATS_STM_H

#include "types.hpp"
#include "la.hpp"
#include "atomistic/fundamental.hpp"
#include "formats/cube.hpp"

namespace formats {
namespace stm {

/**
 * STS in 2d
 * 
 */
class STS2d : public formats::Cube {
    private:
        //std::vector< formats::WfnCube > &levels;
        types::Real broadening;

    public:
        STS2d(
                std::list< formats::WfnCube > &cubes,
                types::Real height,
                types::Real eMin,
                types::Real eMax,
                types::Real deltaE,
                types::Real broadening);
       void addLevel(const std::vector<types::Real> &plane,
               types::Real energy);
       types::Real getEMin();
       types::Real getEMax();
       types::Real getDeltaE();
};
    
    
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
}
#endif
