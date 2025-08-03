#include <iostream>
#include <random>

#include "LinAlg/LinAlg.hpp"

using namespace Planck;

int main() {
    matrix A = {{1,2}, {5,6}};
    matrix B = {{1,3,2,4}, {3,5,4,6}};
    vector b = {1,2,3,4};

    // complex C = dot(B, B);

    // vector d = hermitian(b);

    matrix d = kron(A, B);
    // std::cout << A.size() << " " << A[0].size() << "\n";
    // std::cout << B.size() << " " << B[0].size() << "\n";
    // matrix d = dot(A, B);

    for(int i = 0; i < d.size(); i++){
        for(int j = 0; j < d[i].size(); j++){
            std::cout << d[i][j] << " ";
        }
        std::cout << "\n";
    }

    // for(int i = 0; i < d.size(); i++) std::cout << b[i] << "\n";


    return 0; 
}