#include "QAOA.hpp"
#include "Circuit/Circuit.hpp"
#include "Hamiltonian/Hamiltonian.hpp"
#include "Random/Random.hpp"

#include <cmath>

namespace Planck {

    QAOA::QAOA(Circuit ansatz, matrix Q, vector c) : VQE(ansatz, QUBOHamiltonian(Q, c)), Q(Q), c(c) {

    }

    // autoamtically construct a valid QAOA circuit --
    // set ansatz param to this circuit. 
    // Dummy circuit
    QAOA::QAOA(int p, matrix Q, vector c) : VQE(Circuit(Q.size()), QUBOHamiltonian(Q, c)), Q(Q), c(c) {

        std::cout << "hi kid" << "\n";

        int numQubits = Q.size();
        
        Circuit ansatz = Circuit(Q.size());

        for(int i = 0; i < numQubits; i++) ansatz.addGate(i, "hadamard");
        for(int k = 0; k < p; k++){
            // add cost layer 

            for(int i = 0; i < numQubits; i++){
                double w = 0.5 * c[i].real();
                for(int j = 0; j < Q[i].size(); j++) w += 0.5 * Q[i][j].real();
                ansatz.addGate(i, "rZ", uniform(0, 2 * M_PI), w);
            }

            for(int i = 0; i < Q.size(); i++){
                for(int j = i + 1; j < Q[i].size(); j++){
                    if(Q[i][j].real() == 0) continue;
                    double w = 0.25 * Q[i][j].real();

                    ansatz.addGate(j, i, "CNOT");
                    ansatz.addGate(j, "rZ", uniform(0, 2 * M_PI), w);
                    ansatz.addGate(j, i, "CNOT");
                }
            }


            // add mixer layer

            for(int i = 0; i < numQubits; i++){
                ansatz.addGate(i, "rX", uniform(0, M_PI), 2);
            }
        }
    
        this->ansatz = ansatz;
    }

    double QAOA::evaluate() {
        double bestCost = -std::numeric_limits<double>::infinity();
        std::vector<int> bestBitstring;
    
        for (int i = 0; i < 1000; i++) {
            std::vector<int> bitStr = ansatz.sampleCircuit();
            vector x = zeros(bitStr.size());
            for (int j = 0; j < x.size(); j++) x[j] = bitStr[j];
    
            double currCost = (dot(dot(x, Q), x) + dot(c, x)).real(); // assumes Max-Cut Hamiltonian
            if (currCost > bestCost) {
                bestCost = currCost;
                bestBitstring = bitStr;
            }
        }
    
        std::cout << "Best bitstring: ";
        for (int i = 0; i < bestBitstring.size(); i++) {
            std::cout << bestBitstring[i] << " ";
        }
        std::cout << "\nBest Max-Cut value: " << bestCost << std::endl;
    
        return bestCost;
    }
    

}