#ifndef VQE_HPP
#define VQE_HPP

#include "Circuit/Circuit.hpp"
#include "Constants/Constants.hpp"

namespace Planck {

    class VQE {
        public: 

            VQE(Circuit ansatz, matrix H);

            double expectationValue();

            void gradientDescent(int epochs, double learningRate, bool debug = false);

            int numParams();

            Circuit ansatz;
            matrix H; 

        private:
            vector calculateGradients();
    };
}

#endif