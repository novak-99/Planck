#include "Constants/Constants.hpp"
#include "Circuit/Circuit.hpp"
#include "Bits/Bits.hpp"
#include <string>

namespace Planck {

    class QPE {
        public:
            QPE(vector state, std::string U, matrix A, int m);
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