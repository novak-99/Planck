#include "QPE.hpp"

namespace Planck {

    QPE::QPE(vector state, std::string U, matrix A, int m) : state(state), U(U), A(A), m(m)
    {
        int n = state.size(); 
        Circuit circuit = Circuit(m + n);

        for(int i = 0; i < m; i++){
            circuit.addGate(i, "hadamard");
        }

        std::vector<int> ancilla; 
        std::vector<int> eigenState; 
        for(int i = 0; i < m; i++) ancilla.push_back(i);
        for(int i = m; i < m + n; i++) eigenState.push_back(i);

        // U gates for HHL.
        // just assume it's for HHL. fix later.

        for(int i = 0; i < m; i++){
            // controlled U.
            circuit.addGate(i, eigenState, "CNOT");
            circuit.addTrotter(A, std::pow(2, i), m);
            circuit.addGate(i, eigenState, "CNOT");
        }

        circuit.addGate(ancilla, "IQFT");

        vector zeroState(std::pow(2, m)); 
        zeroState[0] = 1; 
        initial = kron(state, zeroState);
    }

    vector QPE::measureState(){
        return convertBitStrToVector(circuit.sampleCircuit(initial));
    }

    double QPE::evaluate(){
        double sample = convertBitStrToDecimal(circuit.sampleCircuit(initial));
        return sample / std::pow(2, m);
    }

    Circuit QPE::toCircuit(){
        return circuit;
    }



}