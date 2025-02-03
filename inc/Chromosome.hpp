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

    Chromosome(const std::vector<int>& genes);

    Chromosome(size_t genes_size);

    Chromosome(const std::vector<int>& first_half, const std::vector<int>& second_half);

    Chromosome(const Chromosome& chromosome);
	
    ~Chromosome() = default;
	
    Chromosome& operator=(const Chromosome& chromosome);
	
    friend std::ostream& operator<<(std::ostream& os, const Chromosome& chromosome);
};


#endif
