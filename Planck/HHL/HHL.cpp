#include "HHL.hpp"
#include "QPE/QPE.hpp"
#include "IQPE/IQPE.hpp"
#include "Measure/Measure.hpp"
#include <cmath>

namespace Planck {
    HHL::HHL(matrix A, vector b, int m) : A(A/A.norm()), b(b/b.norm()), m(m) {

        int n = b.size();

        Circuit circuit = new Circuit(1);
        QPE qpe(b, "HLL", A, m);

        vector eigenValueState = qpe.measureState();

        circuit = concatenate(circuit, qpe.toCircuit());

        // do controlled rotation 
        
        double lambda = convertBitStrToDecimal(measures(eigenValueState));
        double theta = 2 * std::asin(1/lambda);

        std::vector<int> phaseGates; 
        for(int i = 1; i < m + 1; i++) phaseGates.push_back(i);

        circuit.addGate(0, "rY", theta/2, false);
        circuit.addGate(phaseGates, 0, "CNOT");

        circuit.addGate(0, "rY", -theta/2, false);
        circuit.addGate(phaseGates, 0, "CNOT");

        IQPE iqpe(b, "HLL", A, m);
        circuit = horzcat(circuit,  iqpe.toCircuit());
    }

    vector HHL::evaluate(){
        while(true){
            vector state = circuit.evaluate();
            std::vector<int> sample = measure(state);
            if(sample[sample.size() - 1] == 1){
                return state;
            }
        }
    }
}