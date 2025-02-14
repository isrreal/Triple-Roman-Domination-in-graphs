#ifndef UTIL_FUNCTIONS_HPP
#define UTIL_FUNCTIONS_HPP

#include <vector>
#include <random>
#include "Graph.hpp"
#include "Chromosome.hpp"

size_t getRandomInt(size_t, size_t); 

float getRandomFloat(float, float);

bool feasible(const Graph&, const std::vector<int>&);

bool activeNeighborhoodIsFeasible(const Graph& graph, const std::vector<int>& solution, const std::vector<int>& active_neighbors);

bool feasible(const Graph&, const std::vector<int>&, size_t);

Chromosome& feasibilityCheck(const Graph& , Chromosome&);

void feasibilityCheck(const Graph& graph, std::vector<int>& solution);

void decreaseLabels(const Graph&, std::vector<int>&);

void decreaseLabel(const Graph&, std::vector<int>&, size_t);

Chromosome& fitness(Chromosome&);

int computeRightLowerBound(const Graph&, int);

int computeRightUpperBound(Graph&, int);

#endif
