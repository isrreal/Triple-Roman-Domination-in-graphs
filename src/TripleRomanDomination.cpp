#include "TripleRomanDomination.hpp"

TripleRomanDomination::~TripleRomanDomination() {
    delete this->geneticAlgorithm;
    delete this->ACO;
}

Graph& TripleRomanDomination::getGraph() {
    return this->graph;
}

std::vector<int> TripleRomanDomination::getSolutionACO() {
	return this->solutionACO;
}

std::vector<int> TripleRomanDomination::getSolutionGeneticAlgorithm() {
	return this->solutionGeneticAlgorithm;
}

size_t TripleRomanDomination::getGeneticAlgorithmBestFitness() {
    return this->geneticAlgorithmBestFitness;
}

size_t TripleRomanDomination::getACOBestFitness() {
    return this->ACOBestFitness;
}

void TripleRomanDomination::setNeighbor2(const Graph& graph, Chromosome& solution) {
	std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_int_distribution<int> gap(0, graph.getOrder() - 1);
    
	for (size_t i = 0; i < graph.getOrder(); ++i) {
	    if (solution.genes[i] == 2) {
	        bool onlyLabel0 = true;
	        for (const auto& neighbor : graph.getAdjacencyList(i)) {	
	            if (solution.genes[neighbor] >= 2) {
	                onlyLabel0 = false;
	                break;
	            }
	        }

	        if (onlyLabel0) {
	            auto neighbors = graph.getAdjacencyList(i);
	            size_t randomNeighbor = gap(seed) % neighbors.size();
	            auto it = std::next(neighbors.begin(), randomNeighbor);
	            solution.genes[*it] = 2;
	        }
	    }
	}
}

void TripleRomanDomination::toggleLabels(const Graph& graph, Chromosome& solution) {
	size_t initLabel = 0;
	
	for (size_t i = 0; i < graph.getOrder(); ++i) {
	    if (solution.genes[i] == 4 || solution.genes[i] == 3) { 
			initLabel = solution.genes[i];
			solution.genes[i] = 0;
				
			if (!feasible(graph, solution.genes)) {
	    		solution.genes[i] = 2;
	    		
	    		if (!feasible(graph, solution.genes)) {
	    			solution.genes[i] = 3;
	        		
	        		if (!feasible(graph, solution.genes)) 
	            		solution.genes[i] = initLabel;        
	    		}
			}
		}
	}
}

/**
 * @brief Runs the genetic algorithm and calculates the triple Roman domination number (Gamma3r).
 * 
 * Executes the genetic algorithm with predefined heuristics and computes the total sum of genes
 * in the resulting chromosome to calculate Gamma3r.
 * 
 * @return size_t The calculated triple Roman domination number.
 */
 
void TripleRomanDomination::runGeneticAlgorithm(short int heuristic, bool hasRVNS) {  
    this->geneticAlgorithmBestFitness = 0;

    Chromosome (*selectedHeuristic)(Graph) = nullptr;
    Chromosome (*selectedHeuristicRVNS)(Graph, Chromosome&) = nullptr;

    if (heuristic == 2) {
        selectedHeuristic = heuristic2;
        selectedHeuristicRVNS = heuristic2RVNS;
    } 
    else if (heuristic == 3) {
        selectedHeuristic = heuristic3;
        selectedHeuristicRVNS = heuristic3RVNS;
    } 
 
    else {
        selectedHeuristic = heuristic1;
        selectedHeuristicRVNS = heuristic1RVNS;
    }

    if (hasRVNS) 
        this->geneticAlgorithm->run1(geneticAlgorithm->getGenerations(), selectedHeuristic, heuristic1RVNS);
    else 
        this->geneticAlgorithm->run2(geneticAlgorithm->getGenerations(), selectedHeuristic);

    solutionGeneticAlgorithm = this->geneticAlgorithm->getBestSolution();
    std::for_each(solutionGeneticAlgorithm.begin(), solutionGeneticAlgorithm.end(), [&](int element) {
        this->geneticAlgorithmBestFitness += element;
    });
}


void TripleRomanDomination::runACO() {
   this->ACOBestFitness = 0;
   
   this->ACO->run();
   
   solutionACO = this->ACO->getBestSolution();

   std::for_each(solutionACO.begin(), solutionACO.end(), [&](int element) {
      this->ACOBestFitness += element;
   });
}

bool TripleRomanDomination::feasible(const Graph& graph, std::vector<int> solution) {
    bool isValid = false;

    for (size_t i = 0; i < solution.size(); ++i) {

        isValid = false;

        if (solution[i] == 0) {                                 
            size_t countNeighbors2 = 0;
            size_t countNeighbors3 = 0;
            
            for (auto& neighbor : graph.getAdjacencyList(i)) {          

                if ((countNeighbors2 == 1 && solution[neighbor] >= 3) ||
                    (countNeighbors2 == 2 && solution[neighbor] >= 2)) {
                    isValid = true;
                    break;
                }

                if (countNeighbors3 == 1 && solution[neighbor] >= 2) {
                     isValid = true;
                     break;
                }

                if (solution[neighbor] == 4) {
                    isValid = true;
                    break;
                }
                
                if (solution[neighbor] == 3) 
                    ++countNeighbors3;
                    
                if (solution[neighbor] == 2) 
                    ++countNeighbors2;
            }
            
            if (!isValid)
                return false;
        } 

        else if (solution[i] == 2) {
            bool hasNeighborAtLeast2 = false;
            for (auto& neighbor : graph.getAdjacencyList(i)) {
                if (solution[neighbor] >= 2) {
                    hasNeighborAtLeast2 = true;
                    break; 
                }
            }

            if (!hasNeighborAtLeast2) 
                return false;           
        }
    }
    
    return true;    
}

std::vector<int> TripleRomanDomination::feasibilityCheck(const Graph& graph, std::vector<int> solution) {  
    bool isValid = false;
                                                                                 
    for (size_t i = 0; i < solution.size(); ++i) {
        if (solution[i] == 0) {                                 
            size_t countNeighbors2 = 0;
            size_t countNeighbors3 = 0;
            
            for (auto& neighbor : graph.getAdjacencyList(i)) {          
                                                                                 
                if ((countNeighbors2 == 1 && solution[neighbor] >= 3) ||
                    (countNeighbors2 == 2 && solution[neighbor] >= 2)) {
                    isValid = true;
                    break;
                }
                                                                                 
                if (countNeighbors3 == 1 && solution[neighbor] >= 2) {
                     isValid = true;
                     break;
                }
                                                                                 
                if (solution[neighbor] == 4) {
                    isValid == true;
                    break;
                }
                
                if (solution[neighbor] == 3) 
                    ++countNeighbors3;
                    
                if (solution[neighbor] == 2) 
                    ++countNeighbors2;
            }
            
            if (!isValid) {
            	if (countNeighbors3 == 0 && countNeighbors2 == 0)
                	solution[i] = 3; 
                else if (countNeighbors3 == 1 || countNeighbors2 == 2) 
                    solution[i] = 2;
            }
        }

        else if (solution[i] == 2) {
            bool hasNeighborAtLeast2 = false;
            for (auto& neighbor: graph.getAdjacencyList(i)) {
                if (solution[neighbor] >= 2) {
                    hasNeighborAtLeast2 = true;
                    break; 
                }
            }
                                                                                 
            if (!hasNeighborAtLeast2) 
                solution[i] = 3;           
        }
    }
    
    return solution;
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
 
Chromosome TripleRomanDomination::heuristic1(Graph graph) {
    Chromosome solution(Chromosome(graph.getOrder()));
    std::vector<int> validVertices;
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    
    size_t graphOrder = graph.getOrder();
    size_t choosenVertex = 0;
	
	Graph temp = graph;
	
    while (graph.getOrder() > 0) {
    
       	validVertices.clear();
    	for (const auto& pair : graph.getAdjacencyList()) 
        	validVertices.push_back(pair.first);

        std::uniform_int_distribution<int> gap(0, validVertices.size() - 1);
        choosenVertex = validVertices[gap(seed)];

        if (graph.getVertexDegree(choosenVertex) == 0)
            solution.genes[choosenVertex] = 3;
        else {
            solution.genes[choosenVertex] = 2;

            for (const auto& neighbor : graph.getAdjacencyList(choosenVertex)) {
                if (solution.genes[neighbor] == -1)
                    solution.genes[neighbor] = 0;
            }
        }

        graph.deleteAdjacencyList(choosenVertex);

        for (size_t i = 0; i < graphOrder; ++i) {
            if (graph.vertexExists(i) && graph.getVertexDegree(i) == 0) {
                solution.genes[i] = 3;
                graph.deleteVertex(i);
            }
        }
    }
    
    feasibilityCheck(temp, solution.genes);

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
 
Chromosome TripleRomanDomination::heuristic2(Graph graph) {
    Chromosome solution(Chromosome(graph.getOrder()));
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::vector<size_t> validVertices;
    
    size_t choosenVertex = 0;
    size_t graphOrder = graph.getOrder();
    Graph temp = graph;
	
    while (graph.getOrder() > 0) {
     	validVertices.clear();
        for (const auto& pair : graph.getAdjacencyList()) 
            validVertices.push_back(pair.first);

        std::uniform_int_distribution<int> gap(0, validVertices.size() - 1);
        
        choosenVertex = validVertices[gap(seed)];
		
        if (graph.getVertexDegree(choosenVertex) == 0)
        	solution.genes[choosenVertex] = 3;
        	
        else {
        	solution.genes[choosenVertex] = 4;

		    for (const auto& it: graph.getAdjacencyList(choosenVertex)) {
		        if (solution.genes[it] == -1)
		            solution.genes[it] = 0;
		    }
		}
		
        graph.deleteAdjacencyList(choosenVertex);

        for (size_t i = 0; i < graphOrder; ++i) {
            if (graph.vertexExists(i) && graph.getVertexDegree(i) == 0) {
                solution.genes[i] = 3;
                graph.deleteVertex(i);
            }
        }   
    }
    
    toggleLabels(temp, solution);

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
 
Chromosome TripleRomanDomination::heuristic3(Graph graph) {
    Chromosome solution(Chromosome(graph.getOrder()));
    std::vector<size_t> sortedVertices(graph.getOrder());
    size_t graphOrder = graph.getOrder();
    
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_int_distribution<int> gap(0, graph.getOrder() - 1);
	
	Graph temp = graph;

    for (size_t i = 0; i < graph.getOrder(); ++i)
        sortedVertices[i] = i;

    std::sort(sortedVertices.begin(), sortedVertices.end(),
        [&](size_t a, size_t b) {
            return graph.getVertexDegree(a) > graph.getVertexDegree(b);
    });

    size_t choosenVertex = 0;

    while ((graph.getOrder() > 0) && (choosenVertex < sortedVertices.size())) {
        if (choosenVertex >= sortedVertices.size()) break;

        while (choosenVertex < sortedVertices.size() && 
                (!graph.vertexExists(sortedVertices[choosenVertex]))) {
            ++choosenVertex;
        }

        if (choosenVertex >= sortedVertices.size()) break;
		
		if (graph.getVertexDegree(choosenVertex) == 0)
        	solution.genes[choosenVertex] = 3;
        	
        else {
        	solution.genes[choosenVertex] = 4;

		    for (const auto& it: graph.getAdjacencyList(choosenVertex)) {
		        if (solution.genes[it] == -1)
		            solution.genes[it] = 0;
		    }
		}

        graph.deleteAdjacencyList(sortedVertices[choosenVertex++]);

        for (size_t i = 0; i < graphOrder; ++i) {
            if (graph.vertexExists(i)) {
                if (graph.getVertexDegree(i) == 0) {
                    solution.genes[i] = 3;
                    graph.deleteVertex(i);
                }
            }
        }
    }

    return solution;
}

Chromosome TripleRomanDomination::heuristic1RVNS(Graph graph, Chromosome& chromosome) {
    Chromosome solution(Chromosome(graph.getOrder()));
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::vector<int> destroyedVertices;
    Graph temp = graph;

    for (size_t i = 0; i < graph.getOrder(); ++i) {
        if (chromosome.genes[i] == -1) 
            destroyedVertices.push_back(i);       
    }

    size_t choosenVertex = 0;

    while (destroyedVertices.size() > 0) {
    
		std::uniform_int_distribution<int> gap(0, destroyedVertices.size() - 1);
		choosenVertex = destroyedVertices[gap(seed)];

		if (graph.getVertexDegree(choosenVertex) == 0)
		    solution.genes[choosenVertex] = 3;
		else {
		    solution.genes[choosenVertex] = 2;

		    for (const auto& neighbor : graph.getAdjacencyList(choosenVertex)) {
		        if (solution.genes[neighbor] == -1) 
		            solution.genes[neighbor] = 0;              
		    }
		}

		destroyedVertices.erase(
		    std::remove(destroyedVertices.begin(), destroyedVertices.end(), choosenVertex), 
		    destroyedVertices.end()
		);

		for (size_t i = destroyedVertices.size(); i > 0; --i) {
		    size_t index = i - 1;
		    if (graph.vertexExists(destroyedVertices[index]) && 
		        graph.getVertexDegree(destroyedVertices[index]) == 0) {
		        solution.genes[destroyedVertices[index]] = 3;
		        destroyedVertices.erase(destroyedVertices.begin() + index);
		    }
		}
	}
	
    feasibilityCheck(temp, solution.genes);

    return solution;
}

Chromosome TripleRomanDomination::heuristic2RVNS(Graph graph, Chromosome& chromosome) {
    Chromosome solution(Chromosome(graph.getOrder()));
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    
    size_t choosenVertex = 0;
    size_t graphOrder = graph.getOrder();
    Graph temp = graph;
    
    std::vector<int> destroyedVertices;

    for (size_t i = 0; i < graph.getOrder(); ++i) {
        if (chromosome.genes[i] == -1) 
            destroyedVertices.push_back(i);       
    }
    
    std::uniform_int_distribution<int> gap(0, destroyedVertices.size() - 1);
	
    while (destroyedVertices.size() > 0) {
        choosenVertex = destroyedVertices[gap(seed)];
        while (!graph.vertexExists(choosenVertex))
            choosenVertex = gap(seed);

        if (graph.getVertexDegree(choosenVertex) == 0)
        	solution.genes[choosenVertex] = 3;
        	
        else {
        	solution.genes[choosenVertex] = 4;

		    for (const auto& it: graph.getAdjacencyList(choosenVertex)) {
		        if (solution.genes[it] == -1)
		            solution.genes[it] = 0;
		    }
		}

        destroyedVertices.erase(
            std::remove(destroyedVertices.begin(), destroyedVertices.end(), choosenVertex), 
            destroyedVertices.end()
        );

       	for (size_t i = destroyedVertices.size() - 1; i >= 0; --i) {
            if (graph.vertexExists(destroyedVertices[i]) && graph.getVertexDegree(destroyedVertices[i]) == 0) {
                solution.genes[destroyedVertices[i]] = 3;
                destroyedVertices.erase(destroyedVertices.begin() + i);
            }
        }  
    }
    
    toggleLabels(temp, solution);

    return solution;    
}

Chromosome TripleRomanDomination::heuristic3RVNS(Graph graph, Chromosome& chromosome) {
	Chromosome solution(Chromosome(graph.getOrder()));
    std::vector<size_t> sortedVertices(graph.getOrder());
    size_t graphOrder = graph.getOrder();
	
	Graph temp = graph;
	
	std::vector<int> destroyedVertices;

    for (size_t i = 0; i < graph.getOrder(); ++i) {
        if (chromosome.genes[i] == -1) 
            destroyedVertices.push_back(i);       
    }
  
    for (size_t i = 0; i < graph.getOrder(); ++i)
        sortedVertices[i] = i;

    std::sort(sortedVertices.begin(), sortedVertices.end(),
        [&](size_t a, size_t b) {
            return graph.getVertexDegree(a) > graph.getVertexDegree(b);
    });

    size_t choosenVertex = 0;

    while ((destroyedVertices.size() > 0) && (choosenVertex < sortedVertices.size())) {
        if (choosenVertex >= sortedVertices.size()) break;

        while (choosenVertex < sortedVertices.size() && 
                (!graph.vertexExists(sortedVertices[choosenVertex]))) {
            ++choosenVertex;
        }

        if (choosenVertex >= sortedVertices.size()) break;
		
		if (graph.getVertexDegree(choosenVertex) == 0)
        	solution.genes[choosenVertex] = 3;
        	
        else {
        	solution.genes[choosenVertex] = 4;

		    for (const auto& it: graph.getAdjacencyList(choosenVertex)) {
		        if (solution.genes[it] == -1)
		            solution.genes[it] = 0;
		    }
		}
		
		destroyedVertices.erase(
				    std::remove(destroyedVertices.begin(), destroyedVertices.end(), destroyedVertices[sortedVertices[choosenVertex++]]), 
				    destroyedVertices.end()
		);

        for (size_t i = destroyedVertices.size() - 1; i >= 0; --i) {
            if (graph.vertexExists(destroyedVertices[sortedVertices[choosenVertex++]]) 
            	&& graph.getVertexDegree(destroyedVertices[sortedVertices[choosenVertex++]]) == 0) {
                solution.genes[destroyedVertices[i]] = 3;
                destroyedVertices.erase(destroyedVertices.begin() + i);
            }
        }
    }

    return solution;
}
