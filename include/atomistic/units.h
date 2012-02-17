/*
 * units.h
 */
#ifndef ATOMISTIC_UNITS_H
#define ATOMISTIC_UNITS_H

#include <gsl/gsl_const_mksa.h>
#include "types.hpp"

namespace atomistic {

    namespace units {

        const types::Real eV = GSL_CONST_MKSA_ELECTRON_VOLT;
        const types::Real Ry = GSL_CONST_MKSA_RYDBERG;
        const types::Real Ha = 2 * Ry;
    }
}

#endif
