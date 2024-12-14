#include "Chromosome.hpp"

Chromosome::Chromosome(std::vector<int> genes) {
	this->genes_size = genes.size();
	this->genes = genes;
	this->fitness  = 0;
}

Chromosome::Chromosome(size_t genes_size) {
    this->genes_size = genes_size;
    this->genes = std::vector<int>(genes_size, -1);
    this->fitness  = 0;
} 

Chromosome::Chromosome(std::vector<int> first_half, std::vector<int> second_half) {
    this->genes_size = first_half.size() + second_half.size(); 
    this->genes = first_half;
    this->genes.insert(this->genes.end(), second_half.begin(), second_half.end());
    this->fitness  = 0;
}

Chromosome::Chromosome(const Chromosome& chromosome) {
    this->genes_size = chromosome.genes_size;
    this->genes = chromosome.genes;
    this->fitness  = chromosome.fitness ;
}

std::ostream& operator<<(std::ostream& os, const Chromosome& chromosome) {
    for (const auto& it: chromosome.genes) {
        os << it << " ";
    }
    return os;
}
