#include <vector>
#include <iostream>
#include <string>

#include "stm.hpp"
#include "types.hpp"

namespace stm {

using namespace types;

bool Simulation::getWfnsForBiasRange{
    std::vector< atomistic::WfnCube > wfns;
    
    for(Uint nSpin = 0; nSpin < spectrum.spins.size(); ++nSpin){
        for(Uint nLevel = 0; nLevel < spectrum.spins[nSpin].levels.size(); ++nLevel){
            Real level = spectrum.spins[nSpin].levels[nLevel];
                if(level >= biasLow && level < biasHigh){
                    wfncubes.push_back[


                    binning[nBin].push_back(.size()
               
               

}
