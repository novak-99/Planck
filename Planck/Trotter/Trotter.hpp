#ifndef TROTTER_HPP
#define TROTTER_HPP

#include "Constants/Constants.hpp"
#include "LinAlg/LinAlg.hpp"
#include "Random/Random.hpp"
#include "Gate/Gate.hpp"

class Planck {
    inline std::vector<std::pair<double, std::vector<std::string>>> pauliDecomposition(matrix A){
        int n = A.size();
        int decompSize = std::pow(4, n);

        std::vector<std::string> pauliGates = {"identity", "pauliX", "pauliY", "pauliZ"};
        std::vector<std::vector<std::string>> allSamples = generateAllSamples(pauliGates, n);
        std::vector<std::pair<double, std::vector<string>>> decomp;

        for(int i = 0; i < allSamples.size(); i++){
            matrix currPauliGate = {{1}};
            std::vector<matrix> currPauliString;
            for(int j = 0; j < allSamples[i].size(); j++){
                currPauliGate = kron(currPauliGate, allSamples[i][j]);
                currPauliString.push_back(allSamples[i][j]);
            }
            double coeff = complex(1)/std::pow(2, n) * trace(dot(currPauliGate, A)).real();
            if(coeff != 0) pauliDecomposition.push_back(coeff, currPauliString);
        }

        return decomp;
    }
}

#endif // TROTTER_HPP