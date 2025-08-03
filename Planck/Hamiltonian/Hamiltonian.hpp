#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include "Constants/Constants.hpp"
#include "LinAlg/LinAlg.hpp"

#include <cmath>

namespace Planck{

    inline matrix QUBOHamiltonian(matrix Q, vector c){

        int numQubits = Q.size(); 

        matrix H = zeros(std::pow(2, numQubits), std::pow(2, numQubits));


        for(int i = 0; i < Q.size(); i++){
            for(int j = i + 1; j < Q[i].size(); j++){
                matrix Zij = createGate("pauliZ", {i, j}, numQubits);

                H += 0.25 * Q[i][j] * Zij;

            }
        }


        for(int i = 0; i < Q.size(); i++){
            matrix Zi = createGate("pauliZ", i, numQubits);
            int sumQ = 0; 
            for(int j = 0; j < Q[i].size(); j++) if(i != j) sumQ += Q[i][j].real();
            H -= 0.5 * (c[i].real() + sumQ) * Zi;
        }

        matrix I = eye(std::pow(2, numQubits));

        for(int i = 0; i < Q.size(); i++){
            for(int j = i + 1; j < Q[i].size(); j++){
                H += Q[i][j] / complex(4) * I;
            }
            H += c[i] / complex(2) * I;
        }

        return H;

    }
}

#endif