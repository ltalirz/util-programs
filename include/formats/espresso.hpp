/**
 * @file espresso.hpp
 */
#ifndef FORMATS_ESPRESSO_HPP
#define FORMATS_ESPRESSO_HPP

#include "types.hpp"
#include "atomistic/fundamental.hpp"

namespace formats {

namespace espresso {

/**
 * A whole spectrum, containing several sets of energy levels
 * (spins, kpoints)
 */
struct Spectrum {
    std::vector<atomistic::EnergyLevels> kpoints;

    Spectrum & operator*=(types::Real factor);
    bool readFromOutput(types::String filename);
    void print() const;
    void shift(types::Real deltaE);
    void setFermiZero();
    atomistic::EnergyLevels mergeKpoints() const;
};


}
}
#endif
