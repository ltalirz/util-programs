/**
 * @file cp2k.hpp
 */
#ifndef FORMATS_CP2K_HPP
#define FORMATS_CP2K_HPP

#include "types.hpp"
#include "atomistic/fundamental.hpp"

namespace formats {

namespace cp2k {

/**
 * A whole spectrum, containing several sets of energy levels
 * (spins, kpoints)
 */
struct Spectrum {
    std::vector<atomistic::EnergyLevels> spins;

    Spectrum & operator*=(types::Real factor);
    bool readFromCp2k(types::String filename);
    void print() const;
    void shift(types::Real deltaE);
    void setFermiZero();
    atomistic::EnergyLevels sumSpins() const;
};


}
}
#endif
