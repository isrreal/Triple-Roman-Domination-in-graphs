#ifndef CHROMOSSOME_HPP
#define CHROMOSSOME_HPP

#include <iostream>
#include <random>
#include <vector>
#include "Graph.hpp"

struct Chromosome {
    size_t genesSize;
    std::vector<int> genes;
    size_t fitnessValue;

    Chromosome() = default;

    Chromosome(std::vector<int> genes);

    Chromosome(size_t genesSize);

    Chromosome(std::vector<int> primeiraMetade, std::vector<int> segundaMetade);

    Chromosome(const Chromosome& chromosome);
	
    Chromosome& operator=(const Chromosome& chromosome) {
        if (this != &chromosome) { 
            genesSize = chromosome.genesSize;
            genes = chromosome.genes;
            fitnessValue = chromosome.fitnessValue;
        }
        
        return *this;
    }

    ~Chromosome() = default;

    friend std::ostream& operator<<(std::ostream& os, const Chromosome& chromosome);
};


#endif
