#include "Chromosome.hpp"

Chromosome::Chromosome(std::vector<int> genes) {
	this->genesSize = genes.size();
	this->genes = genes;
	this->indexRemove = 0;
	this->fitnessValue = 0;
}

Chromosome::Chromosome(size_t genesSize) {
    this->genesSize = genesSize;
    this->genes = std::vector<int>(genesSize, -1);
    this->indexRemove = 0;
    this->fitnessValue = 0;
} 

Chromosome::Chromosome(std::vector<int> firstHalf, std::vector<int> secondHalf) {
    this->genesSize = firstHalf.size() + secondHalf.size(); 
    this->genes = firstHalf;
    this->genes.insert(this->genes.end(), secondHalf.begin(), secondHalf.end());
    this->indexRemove = 0;	
    this->fitnessValue = 0;
}

Chromosome::Chromosome(const Chromosome& chromosome) {
    this->genesSize = chromosome.genesSize;
    this->genes = chromosome.genes;
    this->indexRemove = chromosome.indexRemove;
    this->fitnessValue = chromosome.fitnessValue;
}

std::ostream& operator<<(std::ostream& os, const Chromosome& chromosome) {
    for (const auto& it: chromosome.genes)
        os << it << " ";
    return os;
}
