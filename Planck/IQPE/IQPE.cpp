#include "IQPE.hpp"

namespace Planck {

    IQPE::IQPE(vector state, std::string U, matrix A, int m) : state(state), U(U), A(A), m(m)
    {
        int n = state.size(); 
        Circuit circuit = Circuit(m + n);

        std::vector<int> ancilla; 
        std::vector<int> eigenState; 
        for(int i = 0; i < m; i++) ancilla.push_back(i);
        for(int i = m; i < m + n; i++) eigenState.push_back(i);


        circuit.addGate(ancilla, "QFT");

        for(int i = 0; i < m; i++){
            // controlled U.
            circuit.addGate(i, eigenState, "CNOT");
            circuit.addTrotter(A, -std::pow(2, i), m);
            circuit.addGate(i, eigenState, "CNOT");
        }

        for(int i = 0; i < m; i++){
            circuit.addGate(i, "hadamard");
        }

        vector zeroState(std::pow(2, m)); 
        zeroState[0] = 1; 
        initial = kron(state, zeroState);
    }

    vector IQPE::measureState(){
        return convertBitStrToVector(circuit.sampleCircuit(initial));
    }

    double IQPE::evaluate(){
        double sample = convertBitStrToDecimal(circuit.sampleCircuit(initial));
        return sample / std::pow(2, m);
    }

    Circuit IQPE::toCircuit(){
        return circuit;
    }

}