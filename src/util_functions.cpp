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
    bool is_valid {false};
    bool has_neighbor_at_least_2 {false};
    
    for (size_t i {0}; i < solution.size(); ++i) {

        is_valid = false;

        if (solution[i] == 0) {                                 
            size_t count_neighbors_2 {0};
            size_t count_neighbors_3 = {0};
            
            for (auto& neighbor : graph.getAdjacencyList(i)) {          

                if ((count_neighbors_2 == 1 && solution[neighbor] >= 3) ||
                    (count_neighbors_2 == 2 && solution[neighbor] >= 2)) {
                    is_valid = true;
                    break;
                }

                if (count_neighbors_3 == 1 && solution[neighbor] >= 2) {
                     is_valid = true;
                     break;
                }

                if (solution[neighbor] == 4) {
                    is_valid = true;
                    break;
                }
                
                if (solution[neighbor] == 3) { ++count_neighbors_3; }
                    
                if (solution[neighbor] == 2) { ++count_neighbors_2; }
            }
            
            if (!is_valid) { return false; }
        } 

        else if (solution[i] == 2) {
            has_neighbor_at_least_2 = false;
            for (auto& neighbor : graph.getAdjacencyList(i)) {
                if (solution[neighbor] >= 2) {
                    has_neighbor_at_least_2 = true;
                    break; 
                }
            }

            if (!has_neighbor_at_least_2) { return false; }     
        }
    }
    
    return true;    
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
                                                           
    for (size_t i {0}; i < chromosome.genes_size; ++i) {
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
    
    return chromosome;
}

void toggleLabels(const Graph& graph, std::vector<int>& solution) {
	size_t init_label {0};
	
	for (size_t i {0}; i < graph.getOrder(); ++i) {
	    if (solution[i] == 4 || solution[i] == 3) { 
			init_label = solution[i];
			solution[i] = 0;
				
			if (!feasible(graph, solution)) {
	    		solution[i] = 2;
	    		
	    		if (!feasible(graph, solution)) {
	    			solution[i] = 3;
	        		
	        		if (!feasible(graph, solution)) {
	            		solution[i] = init_label;
	            	}        
	    		}
			}
		}
	}
}

Chromosome& fitness(Chromosome& chromosome) {
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
