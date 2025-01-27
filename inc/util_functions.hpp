#ifndef UTIL_FUNCTIONS_HPP
#define UTIL_FUNCTIONS_HPP

#include <vector>
#include <random>
#include "Graph.hpp"
#include "Chromosome.hpp"

size_t getRandomInt(size_t, size_t); 

float getRandomFloat(float, float);

bool feasible(const Graph&, const std::vector<int>&);

Chromosome& feasibilityCheck(const Graph& , Chromosome&);

void toggleLabels(const Graph&, std::vector<int>&);

Chromosome& fitness(Chromosome&);

int computeRightLowerBound(const Graph&, int);

int computeRightUpperBound(Graph&, int);

#endif
