// #ifndef GATE_HPP
// #define GATE_HPP

// #include <iostream>

// #include "Gate.hpp"
// #include "Constants/Constants.hpp"
// #include "LinAlg/LinAlg.hpp"

// namespace Planck {

//     matrix rX(double theta){
//         matrix rX = {
//             {complex(std::cos(theta/2), 0), complex(0, -std::sin(theta/2))},
//             {complex(0, -std::sin(theta/2)), complex(std::cos(theta/2), 0)}
//         };

//         return rX;
//     }

//     matrix rY(double theta){
//         matrix rY = {
//             {complex(std::cos(theta/2), 0), complex(-std::sin(theta/2), 0)},
//             {complex(std::sin(theta/2), 0), complex(std::cos(theta/2), 0)}
//         };
//         return rY;
//     }

//     matrix rZ(double theta){
//         matrix rZ = {
//             {std::exp(complex(0, -theta/2)), 0},
//             {0, std::exp(complex(0, theta/2))}
//         };
//         return rZ;
//     }

//     matrix rX(double theta, double alpha){
//         return rX(theta * alpha);
//     }

//     matrix rY(double theta, double alpha){
//         return rY(theta * alpha);
//     }

//     matrix rZ(double theta, double alpha){
//         return rZ(theta * alpha);
//     }

//     matrix gateManifest(std::pair<std::string, std::vector<double>> gate){
//         if(gate.first == "identity") return identity;
//         else if(gate.first == "pauliX") return pauliX;
//         else if(gate.first == "pauliY") return pauliY;
//         else if(gate.first == "pauliZ") return pauliZ;
//         else if(gate.first == "hadamard") return hadamard;
//         else if(gate.first == "rX") 
//             return gate.second.size() == 1 ? rX(gate.second[0]) : rX(gate.second[0], gate.second[1]);
//         else if(gate.first == "rY") 
//             return gate.second.size() == 1 ? rY(gate.second[0]) : rY(gate.second[0], gate.second[1]);
//         else if(gate.first == "rZ") 
//             return gate.second.size() == 1 ? rZ(gate.second[0]) : rZ(gate.second[0], gate.second[1]);
//     }

//     bool hasCNOT(std::vector<std::pair<std::string, std::vector<double>>> layerGates){
//         for(int i = 0; i < layerGates.size(); i++) {
//             if(layerGates[i].first == "controlCNOT" || layerGates[i].first == "targetCNOT"){
//                 return true; 
//             }
//         }
//         return false;

//     }

//     matrix handleUnary(std::vector<std::pair<std::string, std::vector<double>>> layerGates){
//         matrix layerGate = {{1}};

//         for(int i = 0; i < layerGates.size(); i++){
//             matrix gate = gateManifest(layerGates[i]);


//             layerGate = kron(layerGate, gate);
//         }

//         return layerGate;

//     }

//     matrix handleCNOT(std::vector<std::pair<std::string, std::vector<double>>> layerGates){
//         int control; 
//         int target; 

//         matrix u1 = {{1}};

//         for(int i = 0; i < layerGates.size(); i++){
//             if(layerGates[i].first == "controlCNOT") u1 = kron(u1, p0);
//             else if(layerGates[i].first == "targetCNOT") u1 = kron(u1, identity);
//             else u1 = kron(u1, gateManifest(layerGates[i]));
//         }

//         matrix u2 = {{1}};

//         for(int i = 0; i < layerGates.size(); i++){
//             if(layerGates[i].first == "controlCNOT") u2 = kron(u2, p1);
//             else if(layerGates[i].first == "targetCNOT") u2 = kron(u2, pauliX);
//             else u2 = kron(u2, gateManifest(layerGates[i]));
//         }

//         return u1 + u2;
//     }

//     vector applyLayer(std::vector<std::pair<std::string, std::vector<double>>> layerGates, vector state){
//         matrix layerGate = hasCNOT(layerGates) ? handleCNOT(layerGates) : handleUnary(layerGates);

//         // for(int i = 0; i < layerGate.size(); i++){
//         //     for(int j = 0; j < layerGate[i].size(); j++){
//         //         std::cout << layerGate[i][j] << ' ';
//         //     }
//         //     std::cout << "\n";
//         // }

//         // for(int i = 0; i < state.size(); i++) std::cout << state[i] << '\n';
//         return dot(layerGate, state);
//     } 
// }

// #endif