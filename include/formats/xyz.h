/*
 * xyz.h
 */
#ifndef FORMATS_XYZ_H
#define FORMATS_XYZ_H

#include "types.hpp"
#include "la.h"
#include "atomistic/fundamental.h"

namespace formats {

/**
 * Stores .xyz file
 */
struct Xyz {
    std::vector<atomistic::Atom> atoms;
};

}

#endif
