#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP

#include <string>
#include <utility>
#include <vector>
#include <map>
#include <stdexcept>

#include "Constants/Constants.hpp"
#include "Gate/Gate.hpp"

namespace Planck {

    class Circuit {

        public: 
            Circuit();

            Circuit(int numQubits);
            
            void addQubit();

            void addGate(int qubit, std::string gate);

            // expand these to all use std::vector<T> param later.
            void addGate(int qubit, std::string gate, double param);

            // for no param. 
            // figure out better system for this.
            void addGate(int qubit, std::string gate, double param, bool noParam);

            void addGate(int qubit, std::string gate, double param, double weight);

            void addGate(int control, int target, std::string gate);

            // multi-controlled.
            void addGate(int control, std::vector<int> target, std::string gate);

            void addGate(std::vector<int> control, int target, std::string gate);

            void addTrotter(matrix A, double t, int m);

            // for subcircuits -- like QFT, IQFT, et. al.
            void addGate(std::vector<int> qubits, std::string gate);

            vector evaluate();

            vector evaluate(vector state);

            int numQubits();

            vector getParams();

            void setParams(vector params);

            std::vector<int> sampleCircuit();

            std::vector<int> sampleCircuit(vector state);

            friend Circuit concatenate(Circuit a, Circuit b);
            friend Circuit horzcat(Circuit a, Circuit b);

        private:
            int longestQubit();

            vector params; 
            std::map<std::pair<int, int>, int> paramManifest;
            std::vector<std::vector<std::pair<std::string, std::vector<double>>>> gates;
    };

}

#endif