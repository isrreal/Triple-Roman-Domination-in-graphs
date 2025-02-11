#include "TripleRomanDomination.hpp"

/**
 * @brief Runs the genetic algorithm and calculates the triple Roman domination number (Gamma3r).
 * 
 * Executes the genetic algorithm with predefined heuristics and computes the total sum of genes
 * in the resulting chromosome to calculate Gamma3r.
 * 
 * @return size_t The calculated triple Roman domination number.
 */
 
void TripleRomanDomination::runGeneticAlgorithm(short int heuristic) {  
    this->genetic_algorithm_best_fitness = 0;

    std::vector<std::function<Chromosome(const Graph&)>> heuristics;
    heuristics.reserve(3);
    
    heuristics.emplace_back(heuristic1);
    heuristics.emplace_back(heuristic2);
    heuristics.emplace_back(heuristic3);

   	genetic_algorithm.run(genetic_algorithm.getGenerations(), heuristics, heuristic);
    
    solution_genetic_algorithm = genetic_algorithm.getBestSolution();
    
    this->genetic_algorithm_best_fitness = std::accumulate(solution_genetic_algorithm.begin(), solution_genetic_algorithm.end(), 0);
}

void TripleRomanDomination::runACO(bool with_RVNS) {
   this->aco_best_fitness = 0;
   
   this->ACO.run(with_RVNS);
   
   solution_aco = this->ACO.getBestSolution();

   this->aco_best_fitness = std::accumulate(solution_aco.begin(), solution_aco.end(), 0);   
}

/**
 * @brief A heuristic function that generates an initial chromosome solution for triple Roman domination.
 * 
 * This heuristic randomly selects vertices and assigns them a value of 3, while updating their neighbors 
 * with a value of 0. The adjacency list of the choosen vertex is then deleted.
 * 
 * @param  graph Object used to create the chromosome.
 * @return Chromosome Object The generated chromosome solution.
 */
 
Chromosome TripleRomanDomination::heuristic1(const Graph& graph) {
    Chromosome solution(Chromosome(graph.getOrder()));
    std::vector<int> valid_vertices;
    size_t chosen_vertex {0};
    valid_vertices.reserve(graph.getOrder());
    Graph temp {graph};

    while (temp.getOrder() > 0) {
        valid_vertices.clear();
        for (const auto& pair : temp.getAdjacencyList()) {
            valid_vertices.push_back(pair.first);
        }
        
        chosen_vertex = valid_vertices[getRandomInt(0, valid_vertices.size() - 1)];

        if (temp.getVertexDegree(chosen_vertex) == 0) {
            solution.genes[chosen_vertex] = 3;
        } 
        
        else {
            solution.genes[chosen_vertex] = 2;
            
            for (const auto& neighbor : temp.getAdjacencyList(chosen_vertex)) {
                if (solution.genes[neighbor] == -1) {
                    solution.genes[neighbor] = 0;
                }
            }
        }

        temp.deleteAdjacencyList(chosen_vertex);
        temp.deleteVertex(chosen_vertex);

        std::vector<size_t> vertices_to_remove;
        
        vertices_to_remove.reserve(temp.getOrder());
        
        for (const auto& [i, _]: temp.getAdjacencyList()) {
            if (temp.getVertexDegree(i) == 0) {
                solution.genes[i] = 3;
                vertices_to_remove.push_back(i);
            }
        }

        for (const auto& vertex : vertices_to_remove) {
            temp.deleteVertex(vertex);
        }
    }

    feasibilityCheck(graph, solution);
    fitness(solution);

    return solution;
}


/**
 * @brief A second heuristic function that generates an initial chromosome solution for triple Roman domination.
 * 
 * This heuristic selects vertices randomly and assigns them values, updating neighbors accordingly and
 * handling vertices with no neighbors by assigning a value of 2.
 * 
 * @param graph Object to the graph used to create the chromosome.
 * @return Chromosome Object The generated chromosome solution.
 */
 
Chromosome TripleRomanDomination::heuristic2(const Graph& graph) {
    Chromosome solution(Chromosome(graph.getOrder()));
    std::vector<size_t> valid_vertices;

    size_t chosen_vertex {0};
    
    valid_vertices.reserve(graph.getOrder());
    
    Graph temp {graph};
	
    while (temp.getOrder() > 0) {
     	valid_vertices.clear();
        for (const auto& [i, _] : temp.getAdjacencyList()) {
            valid_vertices.push_back(i);
        }
       
        chosen_vertex = valid_vertices[getRandomInt(0, valid_vertices.size() - 1)];
		
        if (temp.getVertexDegree(chosen_vertex) == 0) {
        	solution.genes[chosen_vertex] = 3;
        }
        	
        else {
        	solution.genes[chosen_vertex] = 4;

		    for (const auto& it: temp.getAdjacencyList(chosen_vertex)) {
		        if (solution.genes[it] == -1) {
		            solution.genes[it] = 0;
		        }
		    }
		}
		
        temp.deleteAdjacencyList(chosen_vertex);
        temp.deleteAdjacencyList(chosen_vertex);
        temp.deleteVertex(chosen_vertex);

       	std::vector<size_t> vertices_to_remove;
        
        vertices_to_remove.reserve(temp.getOrder());
        
        for (const auto& [i, _]: temp.getAdjacencyList()) {
            if (temp.getVertexDegree(i) == 0) {
                solution.genes[i] = 3;
                vertices_to_remove.push_back(i);
            }
        }

        for (const auto& vertex : vertices_to_remove) {
            temp.deleteVertex(vertex);
        }
    }
    
    toggleLabels(graph, solution.genes);
    
    fitness(solution);
	
    return solution;
}


/**
 * @brief A third heuristic function that generates an initial chromosome solution for triple Roman domination.
 * 
 * This heuristic sorts vertices by degree in descending order, then selects vertices with the highest degree
 * for assignment while updating their neighbors.
 * 
 * @param graph Object to the graph used to create the chromosome.
 * @return Chromosome Object The generated chromosome solution.
 */
 
Chromosome TripleRomanDomination::heuristic3(const Graph& graph) {
    Chromosome solution(Chromosome(graph.getOrder()));
    
    std::vector<size_t> sorted_vertices;
    
    sorted_vertices.reserve(graph.getOrder());	
    
    Graph temp {graph};
    
    size_t chosen_vertex {0};
    
   	for (const auto& [i, _] : graph.getAdjacencyList()) {
        sorted_vertices.emplace_back(i);
    }

    std::sort(sorted_vertices.begin(), sorted_vertices.end(),
        [&](size_t a, size_t b) {
            return graph.getVertexDegree(a) < graph.getVertexDegree(b);
    });

    while ((temp.getOrder() > 0) && (chosen_vertex < sorted_vertices.size())) {

        while (chosen_vertex < sorted_vertices.size() && 
            	(!temp.vertexExists(sorted_vertices[chosen_vertex]))) {
            ++chosen_vertex;
        }

        if (chosen_vertex >= sorted_vertices.size()) { break; };

	    if (temp.getVertexDegree(sorted_vertices[chosen_vertex]) == 0) {
	        solution.genes[sorted_vertices[chosen_vertex]] = 3;
	    }
	    
	    else {
	        solution.genes[sorted_vertices[chosen_vertex]] = 4;

	        for (const auto& it : temp.getAdjacencyList(sorted_vertices[chosen_vertex])) {
	            if (solution.genes[it] == -1) {
                	solution.genes[it] = 0;
                }
	        }
	    }

	    temp.deleteAdjacencyList(sorted_vertices[chosen_vertex]);
        temp.deleteVertex(sorted_vertices[chosen_vertex++]); 

	    std::vector<size_t> vertices_to_remove;
        
        vertices_to_remove.reserve(temp.getOrder());
        
        for (const auto& [i, _]: temp.getAdjacencyList()) {
            if (temp.getVertexDegree(i) == 0) {
                solution.genes[i] = 3;
                vertices_to_remove.push_back(i);
            }
        }

        for (const auto& vertex : vertices_to_remove) {
            temp.deleteVertex(vertex);
        }
	}
    
    toggleLabels(graph, solution.genes);
	
    fitness(solution);
	
    return solution;
}

Graph& TripleRomanDomination::getGraph() {
    return this->graph;
}

std::vector<int> TripleRomanDomination::getSolutionACO() {
	return this->solution_aco;
}

std::vector<int> TripleRomanDomination::getSolutionGeneticAlgorithm() {
	return this->solution_genetic_algorithm;
}

size_t TripleRomanDomination::getGeneticAlgorithmBestFitness() {
    return this->genetic_algorithm_best_fitness;
}

size_t TripleRomanDomination::getACOBestFitness() {
    return this->aco_best_fitness;
}
