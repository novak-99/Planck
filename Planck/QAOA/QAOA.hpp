#include "VQE/VQE.hpp"

namespace Planck {
    class QAOA : public VQE {
        public:
            QAOA(Circuit ansatz, matrix Q, vector c); 

            // autoamtically construct a valid QAOA circuit --
            // set ansatz param to this circuit. 
            QAOA(int p, matrix Q, vector c);

            double evaluate();

        private:
            matrix Q; 
            vector c; 

    };
}