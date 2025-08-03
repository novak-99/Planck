#ifndef MEASURE_HPP
#define MEASURE_HPP

#include "LinAlg/LinAlg.hpp"
#include "Random/Random.hpp"
#include "Bits/Bits.hpp"
#include <string>

namespace Planck {

    inline std::vector<int> measure(vector state){
        int numStrings = state.size();
        std::vector<double> amplitudes(numStrings);

        for(int i = 0; i < numStrings; i++){
            amplitudes[i] = std::norm(state[i]);
        }

        return getBitStr(sample(amplitudes), std::log2(numStrings));


    }

}

#endif