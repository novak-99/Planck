#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <vector>
#include <complex>

namespace Planck {
    typedef std::complex<double> complex;
    typedef std::vector<std::vector<complex>> matrix;
    typedef std::vector<complex> vector;

    extern const matrix pauliX;
    extern const matrix pauliY;
    extern const matrix pauliZ;
    extern const matrix hadamard;
    extern const matrix CNOT;
    extern const matrix identity;
    extern const matrix p0;
    extern const matrix p1;
    extern const matrix S; 
    extern const matrix SInv;
}

#endif // CONSTANTS_HPP
