#ifndef BITS_HPP
#define BITS_HPP

#include "Constants/Constants.hpp"
#include "LinAlg/LinAlg.hpp"
#include <vector>

namespace Planck {
    inline std::vector<int> getBitStr(int idx, int bitLength) {
        std::vector<int> bitStr;

        for (int i = bitLength - 1; i >= 0; i--) {
            bitStr.push_back((idx >> i) & 1);
        }

        return bitStr;
    }

    inline vector convertBitStrToVector(std::vector<int> bitStr){
        vector a(bitStr.size());
        for(int i = 0; i < a.size(); i++) a[i] = bitStr[i];
        return a;
    }

    inline int convertBitStrToDecimal(std::vector<int> bitStr){
        int num = 0; 
        int n = bitStr.size();
        for(int i = 0; i < n; i++){
            num += std::pow(2, (n - 1) - bitStr[i]) * bitStr[i];
        }
        return num;
    }
}


#endif // BITS_HPP