
#include "util_functions.hpp"  

std::random_device rd;
std::mt19937 seed(rd());

size_t getRandomInt(size_t start, size_t end) {
    std::uniform_int_distribution<size_t> gap(start, end);
    return gap(seed);
}

float getRandomFloat(float start, float end) {
    std::uniform_real_distribution<float> gap(start, end);
    return gap(seed);
}

bool feasible(const Graph& graph, const std::vector<int>& solution) {
	for (size_t i {0}; i < solution.size(); ++i) {
		if (!feasible(graph, solution, i)) {
			return false;
		}	
	}
	
	return true;
}

bool activeNeighborhoodIsFeasible(const Graph& graph, const std::vector<int>& solution, const std::vector<int>& active_neighbors) {
	for (const auto& vertex: active_neighbors) {
		if (!feasible(graph, solution, vertex)) {
			return false;
		}	
	}
	
	return true;
}

/**
 * @brief Checks if a vertex label in the solution is feasible based on its neighbors.
 * @details This function checks the feasibility of a vertex label in the solution according to the rules 
 * for Triple Roman Domination. The feasibility criteria depend on the current label of the vertex and 
 * the labels of its neighbors. The function works as follows:
 * 
 * - If the vertex has label `0`:
 *   - The vertex is considered valid if it has at least one neighbor with label `4`, or if it has neighbors 
 *     with labels `>=3` or `>=2` in specific patterns.
 * - If the vertex has label `2`:
 *   - The vertex is valid if it has at least one neighbor with label `>=2`.
 * 
 * The function returns `true` if the vertex label is feasible, and `false` otherwise.
 * 
 * @param graph The graph used to check the adjacency of the vertex.
 * @param solution The current solution where vertex labels are stored.
 * @param vertex The index of the vertex to be checked.
 * @return `true` if the vertex label is feasible, `false` otherwise.
 */
 
bool feasible(const Graph& graph, const std::vector<int>& solution, size_t vertex) {
	size_t active {0};
	size_t sum_weight = solution[vertex];	
	for (const auto& it: graph.getAdjacencyList(vertex)) {
		if (solution[it] > 0) {
			++active;
			sum_weight += solution[it];
		}
	}

	if (sum_weight < 3 + active) {
		return false;
	}
	
	else { 
		return true;
	}	
}


/**
 * @brief Checks the feasibility of a chromosome and adjusts if it isn't feasible.
 * 
 * Adjusts genes based on the adjacency list of the graph, ensuring that constraints of Triple Roman Domination are met.
 * 
 * @param chromosome A reference to the chromosome to be checked.
 * @return Chromosome The adjusted chromosome.
 */

Chromosome& feasibilityCheck(const Graph& graph, Chromosome& chromosome) {  
    bool is_valid {false};
    bool has_neighbor_at_least_2 {false};
                                                           
    for (size_t i {0}; i < chromosome.genes.size(); ++i) {		
    	is_valid = false;
 		
        if (chromosome.genes[i] == 0) {    
            size_t count_neighbors_2 {0};
            size_t count_neighbors_3 {0};
            
            for (auto& neighbor : graph.getAdjacencyList(i)) {          
                if (chromosome.genes[neighbor] == 4) {
                 	is_valid = true;
                    break;
                }
                                                         	                        
                if ((count_neighbors_2 == 1 && chromosome.genes[neighbor] >= 3) ||
                    (count_neighbors_2 == 2 && chromosome.genes[neighbor] >= 2)) {
                    is_valid = true;
                    break;
                }
                                                                                 
                if (count_neighbors_3 == 1 && chromosome.genes[neighbor] >= 2) {
                 	is_valid = true;
                 	break;
                }
                                                                                 
                if (chromosome.genes[neighbor] == 3) {
                    ++count_neighbors_3;
                }
                    
                if (chromosome.genes[neighbor] == 2) {
                    ++count_neighbors_2;
                }
            }
            
            if (!is_valid) {
            	if (count_neighbors_3 == 0) {
            		if (count_neighbors_2 == 0) {
            			chromosome.genes[i] = 3;
            		}
            			
            		else if (count_neighbors_2 > 0) {
            			chromosome.genes[i] = 2;
            		}
            	} 
            	
            	else if (count_neighbors_3 == 1) {
            		chromosome.genes[i] = 2;
            	}
            }
        }

        else if (chromosome.genes[i] == 2) {
            has_neighbor_at_least_2 = false;
            for (auto& neighbor : graph.getAdjacencyList(i)) {
                if (chromosome.genes[neighbor] >= 2) {
                    has_neighbor_at_least_2 = true;
                    break; 
                }
            }
                                                                                 
            if (!has_neighbor_at_least_2) {
                chromosome.genes[i] = 3;
            }           
        }
    }
    
    fitness(chromosome);
    
    return chromosome;
}

void feasibilityCheck(const Graph& graph, std::vector<int>& solution) {  
    bool is_valid {false};
    bool has_neighbor_at_least_2 {false};
                                                           
    for (size_t i {0}; i < solution.size(); ++i) {		
    	is_valid = false;
 		
        if (solution[i] == 0) {    
            size_t count_neighbors_2 {0};
            size_t count_neighbors_3 {0};
            
            for (auto& neighbor : graph.getAdjacencyList(i)) {          
                if (solution[neighbor] == 4) {
                 	is_valid = true;
                    break;
                }
                                                         	                        
                if ((count_neighbors_2 == 1 && solution[neighbor] >= 3) ||
                    (count_neighbors_2 == 2 && solution[neighbor] >= 2)) {
                    is_valid = true;
                    break;
                }
                                                                                 
                if (count_neighbors_3 == 1 && solution[neighbor] >= 2) {
                 	is_valid = true;
                 	break;
                }
                                                                                 
                if (solution[neighbor] == 3) {
                    ++count_neighbors_3;
                }
                    
                if (solution[neighbor] == 2) {
                    ++count_neighbors_2;
                }
            }
            
            if (!is_valid) {
            	if (count_neighbors_3 == 0) {
            		if (count_neighbors_2 == 0) {
            			solution[i] = 3;
            		}
            			
            		else if (count_neighbors_2 > 0) {
            			solution[i] = 2;
            		}
            	} 
            	
            	else if (count_neighbors_3 == 1) {
            		solution[i] = 2;
            	}
            }
        }

        else if (solution[i] == 2) {
            has_neighbor_at_least_2 = false;
            for (auto& neighbor : graph.getAdjacencyList(i)) {
                if (solution[neighbor] >= 2) {
                    has_neighbor_at_least_2 = true;
                    break; 
                }
            }
                                                                                 
            if (!has_neighbor_at_least_2) {
                solution[i] = 3;
            }           
        }
    }
}

void decreaseLabels(const Graph& graph, std::vector<int>& solution) {
	for (const auto& vertex : solution) {
	    decreaseLabel(graph, solution, vertex);
	}
}

/**
 * @brief Attempts to change the label of a vertex in the solution while ensuring feasibility.
 * 
 * This function modifies the label of a specific vertex to 0 and checks if the solution 
 * remains feasible. If the solution becomes infeasible, it tries setting the vertex label 
 * to 2, and if still infeasible, to 3. If none of these labels result in a feasible solution, 
 * it reverts the vertex label to its initial value.
 * 
 * @param graph The graph that contains the vertex whose label is being modified.
 * @param solution The vector representing the current labeling of the solution, where each element 
 *                 corresponds to the label of a vertex in the graph.
 * @param vertex The index of the vertex whose label is being decreased.
 * 
 * @note This function modifies the solution vector in place, and ensures that the solution remains 
 *       feasible after modifying the label of the vertex. If no feasible configuration is found, 
 *       the vertex's label is reverted to its initial value.
 */

void decreaseLabel(const Graph& graph, std::vector<int>& solution, size_t vertex) {
	size_t init_label = solution[vertex];
	std::vector<int> active_neighborhood;
	
	
	for (const auto& it: graph.getAdjacencyList(vertex)) {
		active_neighborhood.push_back(it);	
	}
	
	if (init_label == 2) {
		solution[vertex] = 0;
			
		if (!feasible(graph, solution, vertex)) { 
    		solution[vertex] = init_label;      
		}
		
	}
	
	else if (init_label == 3) {
		solution[vertex] = 0;
		
		if (!feasible(graph, solution, vertex)) {
			solution[vertex] = 2;
			if (!feasible(graph, solution, vertex)) {
				solution[vertex] = init_label;      
			}
		}
	}
		
	else if (init_label == 4) {
		solution[vertex] = 0;
		
		if (!feasible(graph, solution, vertex)) {
			solution[vertex] = 2;

			if (!feasible(graph, solution, vertex)) {
				solution[vertex] = 3;
				
				if (!feasible(graph, solution, vertex)) {
		    		solution[vertex] = init_label;      
				}
			}
		}
	}
	
	if (!activeNeighborhoodIsFeasible(graph, solution, active_neighborhood)) {
		solution[vertex] = init_label;      	
	} 
}

Chromosome& fitness(Chromosome& chromosome) {
	chromosome.fitness = 0;
	for (auto& gene: chromosome.genes) {
		chromosome.fitness += gene;
	}
        		
 	return chromosome; 
}

int computeRightLowerBound(const Graph& graph, int lowerBound) {
    lowerBound = -1; 
    if (graph.getMaxDegree() >= 3 && graph.getOrder() >= 2) {
        lowerBound = std::ceil(static_cast<size_t>(4.0 * graph.getOrder() / (graph.getMaxDegree() + 1.0)));
	}
	
    return lowerBound;
}

int computeRightUpperBound(Graph& graph, int upperBound) {
    upperBound = -1;
    auto components { graph.connectedComponents() };
    
    if(components.size() == 1 && graph.getMinDegree() >= 2) {
        upperBound = std::floor(3.0 * graph.getOrder() / 2.0);
    }
    else if (components.size() == 1 && graph.getOrder() >= 3) {
        upperBound = std::floor(7.0 * graph.getOrder() / 4.0);
    }
    else if (components.size() > 1 && graph.getOrder() >= 3) {
        upperBound = 0;
        for (auto& par: components) {
            if(par.first == 1) {
                upperBound += 3;
            }
            else if(par.first == 2) {
                upperBound += 4;
            }            
            else if(par.first >= 3 && par.second >= 2) {
                upperBound += std::floor(static_cast<size_t>(3.0 * par.first / 2.0));
            }           
            else if(par.first >= 3) {
                upperBound += std::floor(static_cast<size_t>(7.0 * par.first / 4.0));
            }
        }
    }

    return upperBound;
}
