/*
 * xyz.h
 */
#ifndef FORMATS_XYZ_H
#define FORMATS_XYZ_H

#include "types.hpp"
#include "la.hpp"
#include "atomistic/fundamental.hpp"

namespace formats {

/**
 * Stores .xyz file
 */
struct Xyz {
    std::vector<atomistic::Atom> atoms;
};

}

#endif
