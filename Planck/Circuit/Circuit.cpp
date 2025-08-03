#include <cmath>
#include <iostream>
#include "Circuit.hpp"

#include "LinAlg/LinAlg.hpp"

#include "Measure/Measure.hpp"

#include "Trotter/Trotter.hpp"

#include <algorithm>

namespace Planck {

    Circuit::Circuit() {

    }

    Circuit::Circuit(int numQubits){
        gates.resize(numQubits);
    }

    void Circuit::addQubit(){ 
        gates.resize(gates.size()+1);
    }

    void Circuit::addGate(int qubit, std::string gate){
        gates[qubit].push_back({gate, {}});
    }

    void Circuit::addGate(int qubit, std::string gate, double param){
        gates[qubit].push_back({gate, {param}});
        
        params.push_back(param);

        paramManifest[{qubit, gates[qubit].size() - 1}] = params.size() - 1;
    }

    void Circuit::addGate(int qubit, std::string gate, double param, bool noParam){
        gates[qubit].push_back({gate, {param}});
    }

    void Circuit::addGate(int qubit, std::string gate, double weight, bool param){
        gates[qubit].push_back({gate, {weight}});
    }

    void Circuit::addGate(int qubit, std::string gate, double param, double weight){
        gates[qubit].push_back({gate, {param, weight}});
        
        params.push_back(param);

        paramManifest[{qubit, gates[qubit].size() - 1}] = params.size() - 1;
    }

    void Circuit::addGate(int control, int target, std::string gate){
        if(gate == "CNOT") {
            int longerQubit = gates[control].size() > gates[target].size() ? control : target;
            int shorterQubit = longerQubit == target ? control : target;
            for(int i = gates[shorterQubit].size(); i < gates[longerQubit].size(); i++) {
                gates[shorterQubit].push_back({"identity", {0}});
            }
            
            gates[control].push_back({"controlCNOT", {0}});
            gates[target].push_back({"targetCNOT", {0}});
        }
    }

    void Circuit::addGate(int control, std::vector<int> target, std::string gate){
        if(gate == "CNOT") {
            // pad affected qubits.
            std::vector<int> qubits = target;
            qubits.push_back(control);
            int longestQubit = qubits[0]; 
            int shortestQubit = qubits[0]; 
            for(int i = 0; i < qubits.size(); i++){
                if(gates[qubits[i]].size() > gates[longestQubit].size()) longestQubit = i; 
                if(gates[qubits[i]].size() < gates[shortestQubit].size()) shortestQubit = i; 
            }


            for(int i = 0; i < qubits.size(); i++){
                for(int j = gates[qubits[i]].size();  j < gates[longestQubit].size(); j++){
                    if(i != longestQubit){
                        gates[qubits[i]].push_back({"identity", {0}});
                    }
                }
            }

            for(int i = 0; i < qubits.size(); i++){
                if(qubits[i] == control) gates[qubits[i]].push_back({"controlCNOT", {0}});
                else gates[qubits[i]].push_back({"targetCNOT", {0}});
            }
        }
    }

    void Circuit::addGate(std::vector<int> control, int target, std::string gate){
        if(gate == "CNOT") {
            // pad affected qubits.
            std::vector<int> qubits = target;
            qubits.push_back(control);
            int longestQubit = qubits[0]; 
            int shortestQubit = qubits[0]; 
            for(int i = 0; i < qubits.size(); i++){
                if(gates[qubits[i]].size() > gates[longerQubit].size()) longerQubit = i; 
                if(gates[qubits[i]].size() < gates[shortestQubit].size()) shortestQubit = i; 
            }


            for(int i = 0; i < qubits.size(); i++){
                for(int j = gates[qubits[i]].size();  j < gates[longestQubit].size(); j++){
                    if(i != longerQubit){
                        gates[qubits[i]].push_back({"identity", {0}});
                    }
                }
            }

            for(int i = 0; i < qubits.size(); i++){
                if(qubits[i] == target) gates[qubits[i]].push_back({"targetCNOT", {0}});
                else gates[qubits[i]].push_back({"controlCNOT", {0}});
            }
        }
    }

    void Circuit::addGate(std::vector<int> qubits, std::string gate){
        if(gate == "QFT"){
            // pad affected qubits.
            // replace this with already implemented multi-cntrl CNOT later.
            int startQubit = *std::min_element(qubits.begin(), qubits.end());
            int endQubit = *std::max_element(qubits.begin(), qubits.end());

            int longestQubit = startQubit;
            int shortestQubit = startQubit;
            for(int i = startQubit; i < endQubit; i++){
                if(gates[i].size() > gates[longerQubit].size()) longerQubit = i; 
                if(gates[i].size() < gates[shortestQubit].size()) shortestQubit = i; 
            }


            for(int i = startQubit; i < endQubit; i++){
                for(int j = gates[i].size();  j < gates[longestQubit].size(); j++){
                    if(i != longerQubit){
                        gates[i].push_back({"identity", {0}});
                    }
                }
            }

            for(int i = startQubit; i < endQubit; i++){
                for(int j = i; j < endQubit; j++){
                    if(j - i == 0) {
                        gates[i].push_back({"hadamard", {0}});
                        for(int k = j; k < endQubit; k++) gates[i].push_back({"identity", {0}});
                        continue; 
                    }
                    double theta = 2 * M_PI / (std::pow(2, j - i + 1));
                    gates[i].push_back({"rZ", {theta/2}});
                    
                    for(int k = j; k < endQubit; k++) gates[k].push_back({"identity", {0}});

                    gates[i].push_back({"targetCNOT", {0}});
                    gates[j].push_back({"controlCNOT", {0}});

                    for(int k = j + 1; k < endQubit; k++) gates[k].push_back({"identity", {0}});

                    gates[i].push_back({"rZ", {-theta/2}});
                    for(int k = j; k < endQubit; k++) gates[k].push_back({"identity", {0}});
                }
            }

            // reverse qubits
            std::reverse(qubits.begin() + startQubit, qubits.end() + endQubit);
        }

        else if(gate == "IQFT"){
            // reverse qubits
            std::reverse(qubits.begin() + startQubit, qubits.end() + endQubit);

            int startQubit = *std::min_element(qubits.begin(), qubit.end());
            int endQubit = *std::max_element(qubits.begin(), qubit.end());

            // padd affected qubits.
            int longestQubit = startQubit;
            int shortestQubit = startQubit;
            for(int i = startQubit; i < endQubit; i++){
                if(gates[i].size() > gates[longestQubit].size()) longestQubit = i; 
                if(gates[i].size() < gates[shortestQubit].size()) shortestQubit = i; 
            }

            for(int i = startQubit; i < endQubit; i++){
                for(int j = gates[i].size();  j < gates[longestQubit].size(); j++){
                    if(i != longestQubit){
                        gates[i].push_back({"identity", {0}});
                    }
                }
            }

            for(int i = endQubit; i >= startQubit; i--){
                for(int j = i; j >= startQubit; j--){
                    if(j - i == 0) {
                        gates[i].push_back({"hadamard", {0}});
                        for(int k = j; k < endQubit; k++) gates[i].push_back({"identity", {0}});
                        continue; 
                    }
                    double theta = 2 * M_PI / (std::pow(2, i - j + 1));
                    gates[i].push_back({"rZ", {theta/2}});
                    
                    for(int k = j; k >= 0; k--) gates[k].push_back({"identity", {0}});

                    gates[i].push_back({"targetCNOT", {0}});
                    gates[j].push_back({"controlCNOT", {0}});

                    for(int k = j - 1; k >= startQubit; k--) gates[k].push_back({"identity", {0}});

                    gates[i].push_back({"rZ", {-theta/2}});
                    for(int k = j; k >= startQubit; k--) gates[k].push_back({"identity", {0}});
                }
            }
        }
    }

    // at some point, add params for t itself in e^iAt. 
    void Circuit::addTrotter(matrix A, double t, int m){
        std::vector<std::pair<double, std::vector<std::string>>> pauliDecomp = pauliDecomposition(A);

        int n = A.size();
        int r = pauliDecomp.size();

        for(int i = 0; i < pauliDecomp.size(); i++){
            std::vector<int> activeQubits; 
            for(int j = m; j < m + n; j++){
                if(pauliDecomp[i].first[j] == "pauliX") { 
                    gates[j].push_back({"hadamard", {0}});
                    activeQubits.push_back(j);
                }
                else if(pauliDecomp[i].first[j] == "pauliY") {
                    gates[j].push_back({"rX", {M_PI/2}});
                    activeQubits.push_back(j); 
                }
                else if(pauliDecomp[i].first[j] == "pauliZ") {
                    gate[j].push_back({"identity", {0}});
                    activeQubits.push_back(j);
                } else gate[j].push_back({"identity", {0}});

            }

            int numActive = activeQubits.size();
            for(int i = 0; i < numActive - 2; i++){
                addGate(activeQubits[i], activeQubits[i+1], "CNOT");
            }

            double theta = t * pauliDecomp[i].second / r; 
            int lastQubit = activeQubits[activeQubits.size() - 1];

            gates[j].push_back({"rZ", {2 * theta}});

            for(int i = numActive - 1; i >= 1; i++){
                addGate(activeQubits[i], activeQubits[i-1], "CNOT");
            }

            for(int j = m; j < m + n; j++){
                if(pauliDecomp[i].first[j] == "pauliX") { 
                    gates[j].push_back({"hadamard", {0}});
                    gates[j].push_back({"identity", {0}});
                }
                else if(pauliDecomp[i].first[j] == "pauliY") {
                    gates[j].push_back({"SInv", {0}});
                    gates[j].push_back({"hadamard", {0}});

                }
                else {
                    gate[j].push_back({"identity", {0}});
                    gates[j].push_back({"identity", {0}});
                }
            }


        }
    }

    // add a padding function later.
    // void Circuit::pad(){

    // }

    vector Circuit::evaluate(){

        vector state = zeros(std::pow(2, numQubits()));
        state[0] = 1;

        for(int i = 0; i < gates[longestQubit()].size(); i++){
            // std::cout << dot(hermitian(state), state)  << "\n";
            std::vector<std::pair<std::string, std::vector<double>>> layerGates;
            for(int j = 0; j < numQubits(); j++){
                i >= gates[j].size() ? layerGates.push_back({"identity", {}}) : layerGates.push_back(gates[j][i]);
            }
            std::reverse(layerGates.begin(), layerGates.end());
            state = applyLayer(layerGates, state);
            // std::cout << state[0] << " " << state[1] << " " << state[2] << " " << state[3] << "\n";
        }
        return state;
    }

    // eventually the above should just be 
    // calling this fn with zero state as the i/p
    vector Circuit::evaluate(vector state){

        for(int i = 0; i < gates[longestQubit()].size(); i++){
            // std::cout << dot(hermitian(state), state)  << "\n";
            std::vector<std::pair<std::string, std::vector<double>>> layerGates;
            for(int j = 0; j < numQubits(); j++){
                i >= gates[j].size() ? layerGates.push_back({"identity", {}}) : layerGates.push_back(gates[j][i]);
            }
            std::reverse(layerGates.begin(), layerGates.end());
            state = applyLayer(layerGates, state);
            // std::cout << state[0] << " " << state[1] << " " << state[2] << " " << state[3] << "\n";
        }
        return state;
    }

    int Circuit::numQubits(){
        return gates.size();
    }

    vector Circuit::getParams() {
        return params;
    }

    void Circuit::setParams(vector params){
        this->params = params;

        for(auto [gateIdx, paramIdx] : paramManifest){
            if(params.size() == 1)
                gates[gateIdx.first][gateIdx.second].second = {this->params[paramIdx].real()};
            else
                gates[gateIdx.first][gateIdx.second].second = {this->params[paramIdx].real(), gates[gateIdx.first][gateIdx.second].second[1]};
        }
    }


    std::vector<int> Circuit::sampleCircuit(){
        return measure(evaluate());
    }

    std::vector<int> Circuit::sampleCircuit(vector state){
        return measure(evaluate(state));
    }

    int Circuit::longestQubit(){
        if(numQubits() == 0) return 0;

        int longest = 0; 

        for(int i = 0; i < numQubits(); i++){
            if(gates[i].size() > gates[longest].size()){
                longest = i; 
            }
        }
        return longest;
    }

   friend Circuit concatenate(Circuit a, Circuit b){
        
        // update the indices for the appended paramManifest.
        std::map<std::pair<int, int>, int> bParamManifestAppended;

        for(auto [gateIdx, paramIdx] : b.paramManifest){
            bParamManifestAppended[{gateIdx.first + n, gateIdx.second}] = b.paramManifest[gateIdx] + a.params.size();
        }

        a.paramManifest.insert(a.paramManifest.end(), bParamManifestAppended.begin(), bParamManifestAppended.end());
        a.params.insert(a.params.end(), b.params.begin(), b.params.end());
        a.gates.insert(a.gates.end(), b.gates.begin(), b.gates.end());

        return a;
   }

   // idx is index of where to insert
   // top is top or bottom, change this to an index later
   friend Circuit horzcat(Circuit a, Circuit b, bool top){

        for(int i = 0; i < a.numQubits() - b.numQubits(); i++) b.addQubit();

        top ? ; : std::reverse(b.gates.begin(), b.gates.end());

        for(int i = 0; i < a.numQubits(); i++) a.gate.insert(a.end(), b.begin(), b.end());

        a.paramManifest.insert(a.paramManifest.begin() + idx, b.paramManifest.begin(), b.paramManifest.end());

        std::map<std::pair<int, int>, int> bParamManifestAppended;

        for(auto [gateIdx, paramIdx] : b.paramManifest){
            if(top)
                bParamManifestAppended[{gateIdx.first, gateIdx.second + a.gates[gateIdx.first].size()}] = b.paramManifest[gateIdx] + a.params.size();
            else
                bParamManifestAppended[{(a.numQubits() - 1) - gateIdx.first, gateIdx.second + a.gates[(a.numQubits() - 1) - gateIdx.first].size()}] = b.paramManifest[gateIdx] + a.params.size();
        }

        a.params.insert(a.params.end(), b.params.begin(), b.params.end());

        return a;
   }
}