#ifndef CHROMOSSOME_HPP
#define CHROMOSSOME_HPP

#include <iostream>
#include <random>
#include <vector>
#include "Graph.hpp"

struct Chromosome {
    size_t genes_size;
    std::vector<int> genes;
    size_t fitness;

    Chromosome() = default;

    Chromosome(std::vector<int> genes);

    Chromosome(size_t genes_size);

    Chromosome(std::vector<int> first_half, std::vector<int> second_half);

    Chromosome(const Chromosome& chromosome);
	
    Chromosome& operator=(const Chromosome& chromosome) {
        if (this != &chromosome) { 
            genes_size = chromosome.genes_size;
            genes = chromosome.genes;
            fitness = chromosome.fitness;
        }
        
        return *this;
    }

    ~Chromosome() = default;

    friend std::ostream& operator<<(std::ostream& os, const Chromosome& chromosome);
};


#endif
