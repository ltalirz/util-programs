/**
 * @file cp.hpp
 */
#ifndef FORMATS_CP_HPP
#define FORMATS_CP_HPP

#include "types.hpp"
#include "atomistic/fundamental.hpp"

namespace formats {

namespace cp {

/**
 * A whole spectrum, containing several sets of energy levels
 * (spins, kpoints)
 */
struct Spectrum {
    std::vector<atomistic::EnergyLevels> spins;

    Spectrum & operator*=(types::Real factor);
    bool readFromCp(types::String filename);
    void print() const;
    void shift(types::Real deltaE);
    void setFermiZero();
    atomistic::EnergyLevels sumSpins() const;
};


}
}
#endif
