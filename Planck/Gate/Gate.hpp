#ifndef GATE_HPP
#define GATE_HPP

#include "Gate.hpp"
#include "Constants/Constants.hpp"
#include "LinAlg/LinAlg.hpp"

#include <cmath>
#include <string>
#include <map>
#include <utility>
#include <iostream>
#include <set>

namespace Planck {

    inline matrix rX(double theta){
        matrix rX = {
            {complex(std::cos(theta/2), 0), complex(0, -std::sin(theta/2))},
            {complex(0, -std::sin(theta/2)), complex(std::cos(theta/2), 0)}
        };

        return rX;
    }

    inline matrix rY(double theta){
        matrix rY = {
            {complex(std::cos(theta/2), 0), complex(-std::sin(theta/2), 0)},
            {complex(std::sin(theta/2), 0), complex(std::cos(theta/2), 0)}
        };
        return rY;
    }

    inline matrix rZ(double theta){
        matrix rZ = {
            {std::exp(complex(0, -theta/2)), 0},
            {0, std::exp(complex(0, theta/2))}
        };
        return rZ;
    }

    inline matrix rX(double theta, double alpha){
        return rX(theta * alpha);
    }

    inline matrix rY(double theta, double alpha){
        return rY(theta * alpha);
    }

    inline matrix rZ(double theta, double alpha){
        return rZ(theta * alpha);
    }

    inline matrix rk(double theta){
        matrix rX = {
           {1, 0},
           {0, std::exp(2 * M_PI * complex(0, 1) / std::pow(2, k))}
        };
        return rk;
    }

    inline matrix gateManifest(std::pair<std::string, std::vector<double>> gate){
        if(gate.first == "identity") return identity;
        else if(gate.first == "pauliX") return pauliX;
        else if(gate.first == "pauliY") return pauliY;
        else if(gate.first == "pauliZ") return pauliZ;
        else if(gate.first == "S") return S;
        else if(gate.first == "SInv") return SInv;
        else if(gate.first == "hadamard") return hadamard;
        else if(gate.first == "rX") 
            return gate.second.size() == 1 ? rX(gate.second[0]) : rX(gate.second[0], gate.second[1]);
        else if(gate.first == "rY") 
            return gate.second.size() == 1 ? rY(gate.second[0]) : rY(gate.second[0], gate.second[1]);
        else if(gate.first == "rZ") 
            return gate.second.size() == 1 ? rZ(gate.second[0]) : rZ(gate.second[0], gate.second[1]);
    }

    inline bool hasCNOT(std::vector<std::pair<std::string, std::vector<double>>> layerGates){
        for(int i = 0; i < layerGates.size(); i++) {
            if(layerGates[i].first == "controlCNOT" || layerGates[i].first == "targetCNOT"){
                return true; 
            }
        }
        return false;

    }

    // inline std::set<int> QFTIndices(std::vector<std::pair<std::string, std::vector<double>>> layerGates){
    //     std::vector<int> indices; 

    //     for(int i = 0; i < layerGates; i++){
    //         if(layerGates[i].first == "QFT"){
    //             indices.push_back(i);
    //         }
    //     }
    //     return indices;
    // }

    inline matrix createGate(std::string gate, int qubit, int numQubits){
        matrix singleton = {{1}};
        for(int i = 0; i < numQubits; i++){
            singleton = (i == qubit) ? kron(singleton, gateManifest({gate, {}})) : kron(singleton, identity);
        }
        return singleton;
    }
    
    inline matrix createGate(std::string gate, std::vector<int> qubits, int numQubits){
        matrix u = {{1}};
        for(int i = 0; i < numQubits; i++){
            bool apply = false; 
            for(int j = 0; j < qubits.size(); j++){
                if(qubits[j] == i) { 
                    apply = true;
                    break;
                }
                u = apply ? kron(u, gateManifest({gate, {}})) : kron(u, identity);
            }
        }
        return u;
    }

    inline matrix handleUnary(std::vector<std::pair<std::string, std::vector<double>>> layerGates){
        matrix layerGate = {{1}};

        for(int i = 0; i < layerGates.size(); i++){
            matrix gate = gateManifest(layerGates[i]);


            layerGate = kron(layerGate, gate);
        }

        return layerGate;

    }

    inline matrix handleCNOT(std::vector<std::pair<std::string, std::vector<double>>> layerGates){
        int control; 
        int target; 

        matrix u1 = {{1}};

        for(int i = 0; i < layerGates.size(); i++){
            if(layerGates[i].first == "controlCNOT") u1 = kron(u1, p0);
            else if(layerGates[i].first == "targetCNOT") u1 = kron(u1, identity);
            else u1 = kron(u1, gateManifest(layerGates[i]));
        }

        matrix u2 = {{1}};

        for(int i = 0; i < layerGates.size(); i++){
            if(layerGates[i].first == "controlCNOT") u2 = kron(u2, p1);
            else if(layerGates[i].first == "targetCNOT") u2 = kron(u2, pauliX);
            else u2 = kron(u2, gateManifest(layerGates[i]));
        }

        return u1 + u2;
    }

    // // add subcircuit.
    // inline std::vector<std::pair<std::string, std::vector<double>>> handleQFTLayer(std::vector<std::pair<std::string, std::vector<double>>> layerGates, std::set<double> qubits, int layer){

    //     matrix qft = {{1}};
    //     int numQubits = state.size();

    //     for(int i = 0; i < numQubits; i++){
    //         for(int j = 0; j < numQubits; j++){

    //         }
    //     }
    // }

    // inline std::vector<std::pair<std::string, std::vector<double>>> handleQFT(std::vector<std::pair<std::string, std::vector<double>>> layerGates){
    //     for(int i = 0; i < layerGates.size(); i++){
    //         for(int j = 0; j < layerGates[i].size(); j++){
    //             if(layerGate[i][j] == "QFT") {
    //                 layerGates = handleQFTLayer(layerGates, layerGates[i][j].second, j);
    //                 break;
    //             }
    //         }
    //     }
    //     return layerGates;
    // }

    inline matrix handleBasicGates(std::vector<std::pair<std::string, std::vector<double>>> layerGates){
        return hasCNOT(layerGates) ? handleCNOT(layerGates) : handleUnary(layerGates);
    }

    inline vector applyLayer(std::vector<std::pair<std::string, std::vector<double>>> layerGates, vector state){

        // sub circuit search

        // std::vector<std::pair<std::string, std::vector<double>>> subLayerGate; 

        // // handle subgates -- otherwise apply basic gate sets like normal.
        // for(int i = 0; i < layerGates.size(); i++){
        //     if(QFTIndices[i].first == "QFT") {
        //         state = dot(handleBasicGates(subLayerGate), state);
        //         state = dot(handleQFT(layerGates[i].second), state);
        //         subLayerGate.clear(); // elif for IQFT, et. al.
        //     } else {
        //         subLayerGate.push_back(layerGates[i]);
        //     }
        // }

        // in case no subcircuits found in the last parts of the circuit
        // use up these final gates.
        // if(handleBasicGates.size() != 0) state = dot(handleBasicGates(subLayerGate), state);

        // layerGates = handleQFT(layerGates);
        // and other subcircuits...

        matrix layerGate = handleBasicGates(layerGates);

        return dot(layerGate, state);
    } 
}

#endif