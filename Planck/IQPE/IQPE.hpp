#ifndef IQPE_HPP
#define IQPE_HPP

#include "Constants/Constants.hpp"
#include "Circuit/Circuit.hpp"
#include "Bits/Bits.hpp"
#include <string>

namespace Planck {

    class IQPE {
        public:
            IQPE(vector state, std::string U, matrix A, int m);
            vector measureState();
            double evaluate();
            Circuit toCircuit();


        private:   
            Circuit circuit; 
            vector state; 
            std::string U;
            matrix A; 
            int m;

            vector initial;
    };
}

#endif // IQPE_HPP