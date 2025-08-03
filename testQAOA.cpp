// {0,72478565 0,72850,0, 0,0,1,40,6 0,60,50,96 0,20,27757 0,340,71,6
//     0,94772976 0,79399593}

#include <iostream>
#include <random>

#include "Circuit/Circuit.hpp"
#include "QAOA/QAOA.hpp"


using namespace Planck;

int main() {

    int p = 2; 
    matrix W = {
        {{0, 1, 1, 1,},
        {1, 0, 1, 0,},
        {1, 1, 0, 1,},
        {1, 0, 1, 0,}}
    };
    matrix Q = -W;
    vector c = {3,2,3,2};

    QAOA optimizer = QAOA(p, Q, c);

    optimizer.gradientDescent(100, 0.1, true);

    std::cout << optimizer.evaluate() << "\n";


    return 0;
}