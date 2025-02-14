#include "Chromosome.hpp"

Chromosome::Chromosome(const std::vector<int>& genes) {
	this->genes = genes;
	this->fitness = 0;
}

Chromosome::Chromosome(size_t genes_size) {
    this->genes = std::vector<int>(genes_size, -1);
    this->fitness = 0;
} 

Chromosome::Chromosome(const std::vector<int>& first_half, const std::vector<int>& second_half) {
    this->genes = first_half;
    this->genes.insert(this->genes.end(), second_half.begin(), second_half.end());
    this->fitness  = 0;
}

Chromosome::Chromosome(const Chromosome& chromosome) {
    this->genes = chromosome.genes;
    this->fitness = chromosome.fitness ;
}

Chromosome& Chromosome::operator=(const Chromosome& chromosome) {
	if (this != &chromosome) { 
	    genes = chromosome.genes;
	    fitness = chromosome.fitness;
	}
	    
	return *this;
}

std::ostream& operator<<(std::ostream& os, const Chromosome& chromosome) {
    for (const auto& it: chromosome.genes) {
        os << it << " ";
    }
    return os;
}
