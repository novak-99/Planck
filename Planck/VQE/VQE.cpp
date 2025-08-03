#include <cmath>
#include <iostream>

#include "VQE.hpp"
#include "LinAlg/LinAlg.hpp"
#include "Utility/Utility.hpp"

#include "Measure/Measure.hpp"

namespace Planck {
    VQE::VQE(Circuit ansatz, matrix H) 
    : ansatz(ansatz), H(H)
    {

    }

    double VQE::expectationValue(){
        vector state = ansatz.evaluate();

        return dot(hermitian(state), dot(H, state)).real();
    }

    void VQE::gradientDescent(int epochs, double learningRate, bool debug){
        for(int i = 0; i < epochs; i++){
            vector params = ansatz.getParams();
            ansatz.evaluate();
            ansatz.setParams(params - learningRate * calculateGradients());

            double expectation =  expectationValue();
            if(debug) std::cout << "Epoch: " << i + 1 << " Expectation Value: " << expectation << "\n";
        }
    }

    vector VQE::calculateGradients(){
        vector gradients = zeros(numParams());
        vector paramsOrig = ansatz.getParams();

        for(int i = 0; i < numParams(); i++){
            vector paramsPlus = paramsOrig; 

            paramsPlus[i] += complex(M_PI / 2, 0); 

            ansatz.setParams(paramsPlus);
            
            double expectationPlus = expectationValue();
            
            vector paramsMinus = paramsOrig; 

            paramsMinus[i] -= complex(M_PI / 2, 0); 
            ansatz.setParams(paramsMinus);

            double expectationMinus = expectationValue();


            gradients[i] = 0.5 * (expectationPlus - expectationMinus);

            ansatz.setParams(paramsOrig);

        }

        return gradients;
    }

    int VQE::numParams(){
        return ansatz.getParams().size();
    }

}

