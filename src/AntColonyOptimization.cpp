#include "AntColonyOptimization.hpp"

void AntColonyOptimization::initializePheromones(std::vector<float>& graph_pheromones) {
    for (auto& vertex_pheromones: graph_pheromones) {
        vertex_pheromones = 0.5;
    }
}

/**
 * @brief Constructs a feasible solution for the Triple Roman Domination problem.
 * 
 * @details This function selects a random vertex and labels it with a value of 4, setting the labels 
 * of all its adjacent vertices to 0. These vertices are then removed from the temporary graph.
 * If any remaining vertices in the temporary graph have a degree of 0, they are labeled with 3 and 
 * also removed. This process continues until there are no vertices left in the temporary graph.
 *
 * @return A vector of integers representing the computed solution for the graph, where each element 
 * indicates the label assigned to each vertex according to the Triple Roman Domination rules.
 */

void AntColonyOptimization::constructSolution(std::vector<int>& solution) {
    Graph temp(this->graph);
    size_t vertex {0};
	
    while (temp.getOrder() > 0) {
        vertex = chooseVertex(temp);  
        solution[vertex] = 4;

        for (const auto& neighbor: temp.getAdjacencyList(vertex)) {
            if (solution[neighbor] == -1) {
                solution[neighbor] = 0;
            }
        }

        temp.deleteAdjacencyList(vertex);
        temp.deleteVertex(vertex);
        
        std::vector<size_t> vertices_to_remove;
        
        vertices_to_remove.reserve(temp.getOrder());
        
        for (const auto& [i, _]: temp.getAdjacencyList()) {
            if (temp.getVertexDegree(i) == 0) {
                solution[i] = 3;
                vertices_to_remove.push_back(i);
            }
        }
		
        for (const auto& vertex : vertices_to_remove) {
            temp.deleteVertex(vertex);
        }
    }  
}

void AntColonyOptimization::constructSolutionRVNS(std::vector<int>& solution) {
    Graph temp(this->graph);
    std::vector<int> destroyed_vertices;
    destroyed_vertices.resize(temp.getOrder());
	
    for (size_t i {0}; i < solution.size(); ++i) {
    	if (solution[i] == -1) {
			solution[i] = 0;
    	}    
    }  
    
	feasibilityCheck(temp, solution);
}

/**
 * @brief Increases the labels of selected vertices in the input solution.
 * 
 * @details This function iteratively selects random vertices from the input graph that are labeled 
 * with 0, 2, or 3, and increases their label to 4. The number of vertices to be updated is determined 
 * by the constant `add_vertices_rate` and the count of vertices with labels 0, 2, or 3. This process 
 * continues until the target number of vertices is updated or there are no more eligible vertices.
 * 
 * @return A vector of integers representing the updated solution, where selected labels have been 
 * increased to 4 based on the specified rate.
 */

void AntColonyOptimization::extendSolution(std::vector<int>& solution) {
    size_t itr {0};
    size_t vertex {0};
    std::vector<int> two_or_zero_labeled_vertices;
    two_or_zero_labeled_vertices.reserve(solution.size());
    
    for (size_t i {0}; i < solution.size(); ++i) {
        if (solution[i] <= 3) {
            two_or_zero_labeled_vertices.push_back(i);
        }
    }

    itr = static_cast<size_t>(add_vertices_rate_extend_solution * two_or_zero_labeled_vertices.size());

    while (itr != 0 && !two_or_zero_labeled_vertices.empty()) {          
        vertex = chooseVertex(two_or_zero_labeled_vertices);         
        solution[vertex] = 4; 
        two_or_zero_labeled_vertices.erase(two_or_zero_labeled_vertices.begin() + vertex);                                                                
        --itr;
    }
}

/**	@brief Tries to reduce the labels of vertices 
*	@details Sorts the vertices of graph on descendent order based on its respective degree,
*			and selects its first vertex (that has higher degree value). After this, 
			tries to reduce the values of vertices in solution that have labels 4 or 3, to 2,
			if the solution are feasible, then the label is updated to 0, otherwise,
			the original label is recovered.
*
* @return A solution that can have its label values reduced
*/

void AntColonyOptimization::reduceSolution(std::vector<int>& solution) {
    Graph temp { this->graph };
    std::vector<int> sorted_vertices;
    sorted_vertices.reserve(temp.getOrder());
    
    size_t chosen_vertex {0};
    
    for (size_t i {0}; i < temp.getOrder(); ++i) {
        sorted_vertices.push_back(i);
    }

    std::sort(sorted_vertices.begin(), sorted_vertices.end(),
        [&](size_t a, size_t b) {
            return temp.getVertexDegree(a) >  
            temp.getVertexDegree(b);                                                                             
        });

    while ((temp.getOrder() > 0) && (chosen_vertex < sorted_vertices.size())) {

        while (chosen_vertex < sorted_vertices.size() && 
            	(!temp.vertexExists(sorted_vertices[chosen_vertex]))) {
            ++chosen_vertex;
        }

        if (chosen_vertex >= sorted_vertices.size()) { break; };

        if (solution[sorted_vertices[chosen_vertex]] == 4 ||
         	solution[sorted_vertices[chosen_vertex]] == 3 ||
         	solution[sorted_vertices[chosen_vertex]] == 2) {

			decreaseLabel(this->graph, solution, sorted_vertices[chosen_vertex]);
		}
		
		temp.deleteAdjacencyList(sorted_vertices[chosen_vertex]);
        temp.deleteVertex(sorted_vertices[chosen_vertex++]);    
    }
}

/**
 * @brief Attempts to improve the final solution computed by the algorithm.
 * 
 * @details This function takes a previously computed solution and applies the following subroutines:
 * - @subroutine destroySolution
 * - @subroutine constructSolutionRVNS
 * - @subroutine extendSolution
 * - @subroutine reduceSolution
 * 
 * After each iteration, the function compares the weight of the newly generated solution with the original one.
 * If the new solution has a lower weight, it updates the best solution to this new value.
 * The process repeats until the maximum number of iterations without improvement is reached.
 * 
 * @param solution A vector of integers representing the previously computed by the algorithm.
 * 
 * @return A vector of integers representing the solution with the least weight found by the algorithm.
 */

void AntColonyOptimization::RVNS(std::vector<int>& solution) {
    size_t current_no_improvement_iteration {0};
    size_t max_iterations { max_rvns_iterations };
    std::vector<int> temp { solution };
    current_rvns_number = 1;
	
    while ((current_no_improvement_iteration < max_rvns_no_improvement_iterations) && (max_iterations > 0)) {
        destroySolution(temp);
        constructSolutionRVNS(temp);
        extendSolution(temp);
        reduceSolution(temp);
		
        if (summation(temp) < summation(solution)) {
            solution = temp;
            current_rvns_number = 1;
            current_no_improvement_iteration = 0;
        }

        else {
            ++current_rvns_number;
            ++current_no_improvement_iteration;

            if (current_rvns_number > max_rvns_functions) {
                current_rvns_number = 1;
            }
        }

        --max_iterations;
    }
}

size_t AntColonyOptimization::chooseVertex(const Graph& temp) {
    
    float number { getRandomFloat(0.0, 1.0) }; 
    
    size_t vertex { getRandomInt(0, this->graph.getOrder() - 1) };
    float value { 0.0 };
    
    while (!temp.vertexExists(vertex)) {
        vertex = getRandomInt(0, this->graph.getOrder() - 1);
    }

    if (selection_vertex_rate_construct_solution < number) {
    	for (const auto& [vertex_iterator, _]: temp.getAdjacencyList()) {
            if (value < (temp.getVertexDegree(vertex_iterator) * graph_pheromones[vertex_iterator])) {
                value = temp.getVertexDegree(vertex_iterator) * graph_pheromones[vertex_iterator];
            	vertex = vertex_iterator;
           	}            																												                                  																																																																											          																																																																					                                                                                                              }
    }

    else {
        vertex = rouletteWheelSelection(temp);
	}
	
    return vertex;
}

/**
 * @brief Selects the index of a vertex that maximizes the objective function or selects one at random.
 * 
 * @details This function either chooses a vertex index that maximizes the product of its degree 
 * and pheromone level (objective function) or randomly selects an index from the 
 * `two_or_zero_labeled_vertices` vector. The choice between these strategies is determined by the 
 * `selection_vertex_rate_extend_solution`, which acts as a threshold probability. If a randomly generated 
 * number exceeds this rate, the function will maximize the objective function; otherwise, a random 
 * index is selected.
 * 
 * @param two_or_zero_labeled_vertices A vector containing indices of vertices labeled as 0, 2, or 3.
 * @return The index within `two_or_zero_labeled_vertices` that either maximizes the objective function
 *         or is randomly selected.
 */

size_t AntColonyOptimization::chooseVertex(const std::vector<int>& two_or_zero_labeled_vertices) {
    size_t chosen_index { getRandomInt(0, two_or_zero_labeled_vertices.size() - 1) };

    if (getRandomFloat(0.0, 1.0) > selection_vertex_rate_extend_solution) {
        float max_objective_value {-1.0f};
        float objective_value {0.0};
        size_t vertex {0};

        for (size_t i {0}; i < two_or_zero_labeled_vertices.size(); ++i) {
            vertex = two_or_zero_labeled_vertices[i];

            objective_value = this->graph.getVertexDegree(vertex) * graph_pheromones[vertex];

            if (objective_value > max_objective_value) {
                max_objective_value = objective_value;
                chosen_index = i;
            }
        }
    }

    else {
        chosen_index = rouletteWheelSelection(two_or_zero_labeled_vertices);
    }

    return chosen_index;
}


size_t AntColonyOptimization::rouletteWheelSelection(const Graph& temp) {
    float total_fitness { 0.0f };
    std::vector<std::pair<size_t, float>> probabilities;
    
    probabilities.reserve(temp.getOrder());

    for (const auto& [i, _]: temp.getAdjacencyList()) {
        total_fitness += temp.getVertexDegree(i) * graph_pheromones[i];
    }

    for (const auto& [i, _]: temp.getAdjacencyList()) {
        probabilities.push_back( {i, static_cast<float>( ((temp.getVertexDegree(i) * graph_pheromones[i]) / total_fitness)) } );
    }

   	float cumulation_sum { 0.0f };
   	
    for (const auto& [vertex, pickRate] : probabilities) {
        cumulation_sum += pickRate;
        if (getRandomFloat(0.0, 1.0) <= cumulation_sum) {
            return vertex;
        }
    }

    return probabilities.back().first;
}



size_t AntColonyOptimization::rouletteWheelSelection(const std::vector<int>& two_or_zero_labeled_vertices) {
    float total_fitness {0.0f};
    std::vector<std::pair<size_t, float>> probabilities;
    
    probabilities.reserve(two_or_zero_labeled_vertices.size());

    for (size_t i {0}; i < two_or_zero_labeled_vertices.size(); ++i)  {
        total_fitness += graph.getVertexDegree(two_or_zero_labeled_vertices[i]) * graph_pheromones[two_or_zero_labeled_vertices[i]];
    }
    
    for (size_t i {0}; i < two_or_zero_labeled_vertices.size(); ++i) {
        probabilities.push_back({i, ((graph.getVertexDegree(two_or_zero_labeled_vertices[i]) * graph_pheromones[two_or_zero_labeled_vertices[i]]) / total_fitness)});
    }
                                             
    float cumulation_sum { 0.0f };
    
    for (const auto& [vertex, pickRate] : probabilities) {
        cumulation_sum += pickRate;
        if (getRandomFloat(0.0, 1.0) <= cumulation_sum) {
            return vertex;
        }
    }
                                             
    return probabilities.back().first;
}


/**
 * @brief Attempts to reduce the labels of vertices in the solution.
 * 
 * @details This function sorts the vertices of the graph in descending order based on their degree 
 * and selects the vertex with the highest degree. For each selected vertex with a label of 4 or 3 
 * in the solution, it attempts to reduce the label to 2 and checks if the solution remains feasible.
 * If the solution remains feasible, the label is updated to 0; otherwise, the original label is restored.
 * This process continues until all vertices have been processed.
 * 
 * @return A vector of integers representing the solution with potentially reduced labels, while 
 * maintaining feasibility according to the problem's constraints.
 */

void AntColonyOptimization::destroySolution(std::vector<int>& solution) {
    float destruction_rate { min_destruction_rate + ((current_rvns_number - 1) *
                ((max_destruction_rate - min_destruction_rate)) 
                / (max_rvns_functions - 1)) };
                
    Graph temp { this->graph };
    size_t itr { static_cast<size_t>(solution.size() * destruction_rate) }; 
    size_t vertex {0};
    
    while (itr != 0 && (temp.getOrder() > 0)) {
       vertex = chooseVertex(temp);
       
       if ((solution[vertex] == 0) || (solution[vertex] == 2) || (solution[vertex] == 3)) {
            solution[vertex] = -1;
       }
       
       else {
            ++itr;
       }
       
       temp.deleteVertex(vertex);
       --itr;
    }
} 

void AntColonyOptimization::updatePheromones(std::vector<int>& current_best_solution, std::vector<int>& best_solution) {   
	short weight_current_best_solution {};
	short weight_best_solution {};
	
	if (convergence_factor < 0.4) {
		weight_current_best_solution = 1;
    	weight_best_solution = 0;
	}                        
	
	else if (convergence_factor >= 0.4 && convergence_factor < 0.6) {
		weight_current_best_solution = static_cast<short>(2 / 3);
    	weight_best_solution = static_cast<short>(1 / 3);
	}        
	
	else if (convergence_factor >= 0.6 && convergence_factor < 0.8) {
		weight_current_best_solution = static_cast<short>(1 / 3);
    	weight_best_solution = static_cast<short>(2 / 3);
	}                   
	
	else {
		weight_current_best_solution = 0;
    	weight_best_solution = 1;
   	}         
    
    float equation {0.0};        
	
    for (size_t vertex {0}; vertex < graph_pheromones.size(); ++vertex) {  
        equation = weight_current_best_solution * delta(current_best_solution, vertex) + 
                      weight_best_solution * delta(best_solution, vertex);
                                                                                                     
        graph_pheromones[vertex] += evaporation_rate * equation;
        
        if (graph_pheromones[vertex] > 0.999) {
        	graph_pheromones[vertex] = 0.999;
        }
        
        if (graph_pheromones[vertex] > 0.001) {
        	graph_pheromones[vertex] = 0.001;
        }
    }
}

float AntColonyOptimization::computeConvergence(const std::vector<float>& graph_pheromones) {
    float max_pheromone { getMaxPheromoneValue(graph_pheromones) }; 
    float min_pheromone { getMinPheromoneValue(graph_pheromones) }; 
    float temp { 0.0 };

    size_t number_of_ants { graph_pheromones.size() };

    for (size_t i {0}; i < number_of_ants; ++i) {
        temp += std::max(max_pheromone - graph_pheromones[i], graph_pheromones[i] - min_pheromone);
    }

    this->convergence_factor = 2 * ((temp / (number_of_ants * (max_pheromone + min_pheromone)))) - 1;

    return this->convergence_factor;
}

size_t AntColonyOptimization::summation(const std::vector<int>& solution) {
    size_t temp {0};
    for (const auto& it: solution) {
        temp += it;
    }
    
    return temp;
}


/* @brief checks if the vertex was choosen
 * The critery of selection vertex is that maximize the function degree(v) * pheromone(v)
 * if 4, so the vertex was choosen
 * 
 *@return the vertex was choosen or not
 *
 */ 
 
bool AntColonyOptimization::delta(const std::vector<int>& solution, size_t vertex) {
    return solution[vertex] == 4 ? true : false; 
}

float AntColonyOptimization::getMinPheromoneValue(const std::vector<float>& graph_pheromones) {
    float value {graph_pheromones[0]};  
                                       
    for (const auto& it: graph_pheromones) {
        if (value > it) {
            value = it;
        }
    }
    
    return value;                   
}

float AntColonyOptimization::getMaxPheromoneValue(const std::vector<float>& graph_pheromones) {
    float value {graph_pheromones[0]};

    for (const auto& it: graph_pheromones) {
        if (value < it) {
            value = it;
        }
    }
    return value;
}

// public methods 

std::vector<int> AntColonyOptimization::getBestSolution() { return this->best_solution; }

void AntColonyOptimization::run(bool with_RVNS) {
    size_t iteration { iterations };
    
    std::vector<int> current_best_solution(graph.getOrder(), 4);
    std::vector<int> best_solution(graph.getOrder(), 4);
    std::vector<int> solution(graph.getOrder(), -1);
	
    initializePheromones(graph_pheromones);
    
    while (iteration > 0) {	  
        for (size_t i {0}; i < number_of_ants; ++i) {
            constructSolution(solution);
			extendSolution(solution);
			reduceSolution(solution);
			
			if (with_RVNS) {
				RVNS(solution);
			}
			
            if (summation(solution) < summation(current_best_solution)) {
                current_best_solution.swap(solution);
            }
        }    
                                           
        if (summation(current_best_solution) < summation(best_solution)) {
              best_solution.swap(current_best_solution);
        }
        
        updatePheromones(current_best_solution, best_solution);
                                                               
        convergence_factor = computeConvergence(graph_pheromones);
                                 
        if (convergence_factor > 0.99) {
            initializePheromones(graph_pheromones);
		}
		
        --iteration; 			
    }
    
    this->best_solution = best_solution;
}
