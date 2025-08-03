#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>
#include <cmath>

namespace Planck {

    inline double uniform(double a, double b){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(a, b); 

        return dist(gen);
    }

    inline double sample(std::vector<double> probabilities){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::discrete_distribution<> dist(probabilities.begin(), probabilities.end()); 

        return dist(gen);
    }

    template<class T>
    inline std::vector<std::vector<T>> generateAllSamples(std::vector<T> dist, int n){
        int d = dist.size();
        int numSamples = std::pow(d, n);

        std::vector<std::vector<T>> samples(numSamples);

        for(int i = 0; i < numSamples; i++){
            int curr = d; 
            for(int j = 0; j < n; j++){
                int idx = curr % d; 
                samples[i].push_back(dist[j]);
                curr /= d; 
            }
        }

        return samples;
    }
}


#endif 