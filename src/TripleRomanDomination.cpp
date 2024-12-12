#include "TripleRomanDomination.hpp"

GeneticAlgorithm* TripleRomanDomination::geneticAlgorithm = nullptr;
AntColonyOptimization* TripleRomanDomination::ACO = nullptr;

TripleRomanDomination::TripleRomanDomination(Graph& graph, size_t populationSize, size_t genesSize, size_t generations, 
    float mutationRate, float elitismRate, float crossOverRate, 
    size_t numberOfAnts, size_t iterations)
    : graph(graph), 
      solutionACO(), 
      solutionGeneticAlgorithm(),
      geneticAlgorithmBestFitness(0), 
      ACOBestFitness(0) {
    	geneticAlgorithm = new GeneticAlgorithm(graph, populationSize, genesSize, generations, mutationRate, elitismRate, crossOverRate);
    	ACO = new AntColonyOptimization(graph, iterations, numberOfAnts);
}



/**
 * @brief Runs the genetic algorithm and calculates the triple Roman domination number (Gamma3r).
 * 
 * Executes the genetic algorithm with predefined heuristics and computes the total sum of genes
 * in the resulting chromosome to calculate Gamma3r.
 * 
 * @return size_t The calculated triple Roman domination number.
 */
 
void TripleRomanDomination::runGeneticAlgorithm(short int heuristic) {  
    this->geneticAlgorithmBestFitness = 0;

    std::vector<std::function<Chromosome(const Graph&)>> heuristics;
    heuristics.reserve(3);
    
    heuristics.emplace_back(heuristic1);
    heuristics.emplace_back(heuristic2);
    heuristics.emplace_back(heuristic3);

    TripleRomanDomination::geneticAlgorithm->run(geneticAlgorithm->getGenerations(), heuristics, heuristic);
    
    solutionGeneticAlgorithm = TripleRomanDomination::geneticAlgorithm->getBestSolution();// heuristic1(graph).genes; 
    
    this->geneticAlgorithmBestFitness = std::accumulate(solutionGeneticAlgorithm.begin(), solutionGeneticAlgorithm.end(), 0);
}

void TripleRomanDomination::runACO() {
   this->ACOBestFitness = 0;
   
   this->ACO->run();
   
   solutionACO = this->ACO->getBestSolution();

   this->ACOBestFitness = std::accumulate(solutionACO.begin(), solutionACO.end(), 0);   
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
    std::vector<int> validVertices;
    size_t choosenVertex {0};
    validVertices.reserve(graph.getOrder());
    Graph temp {graph};

    while (temp.getOrder() > 0) {
        validVertices.clear();
        for (const auto& pair : temp.getAdjacencyList()) {
            validVertices.push_back(pair.first);
        }
        
        choosenVertex = validVertices[GeneticAlgorithm::getRandomInt(0, validVertices.size() - 1)];

        if (temp.getVertexDegree(choosenVertex) == 0) {
            solution.genes[choosenVertex] = 3;
        } 
        
        else {
            solution.genes[choosenVertex] = 2;
            
            for (const auto& neighbor : temp.getAdjacencyList(choosenVertex)) {
                if (solution.genes[neighbor] == -1) {
                    solution.genes[neighbor] = 0;
                }
            }
        }

        temp.deleteAdjacencyList(choosenVertex);

        std::vector<size_t> verticesToRemove;
        
        verticesToRemove.reserve(temp.getOrder());
        
        for (const auto& [i, _]: temp.getAdjacencyList()) {
            if (temp.getVertexDegree(i) == 0) {
                solution.genes[i] = 3;
                verticesToRemove.push_back(i);
            }
        }

        for (const auto& vertex : verticesToRemove) {
            temp.deleteVertex(vertex);
        }
    }

    GeneticAlgorithm::feasibilityCheck(graph, solution);
    GeneticAlgorithm::fitness(solution);

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
    std::vector<size_t> validVertices;

    size_t choosenVertex {0};
    
    validVertices.reserve(graph.getOrder());
    
    Graph temp {graph};
	
    while (temp.getOrder() > 0) {
     	validVertices.clear();
        for (const auto& [i, j] : temp.getAdjacencyList()) {
            validVertices.push_back(i);
        }
       
        choosenVertex = validVertices[GeneticAlgorithm::getRandomInt(0, validVertices.size() - 1)];
		
        if (temp.getVertexDegree(choosenVertex) == 0) {
        	solution.genes[choosenVertex] = 3;
        }
        	
        else {
        	solution.genes[choosenVertex] = 4;

		    for (const auto& it: temp.getAdjacencyList(choosenVertex)) {
		        if (solution.genes[it] == -1) {
		            solution.genes[it] = 0;
		        }
		    }
		}
		
        temp.deleteAdjacencyList(choosenVertex);

       	std::vector<size_t> verticesToRemove;
        
        verticesToRemove.reserve(temp.getOrder());
        
        for (const auto& [i, _]: temp.getAdjacencyList()) {
            if (temp.getVertexDegree(i) == 0) {
                solution.genes[i] = 3;
                verticesToRemove.push_back(i);
            }
        }

        for (const auto& vertex : verticesToRemove) {
            temp.deleteVertex(vertex);
        }
    }
    
    AntColonyOptimization::toggleLabels(graph, solution.genes);
    
    GeneticAlgorithm::fitness(solution);
	
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
    
    std::vector<size_t> sortedVertices;
    
    sortedVertices.reserve(graph.getOrder());	
    
    Graph temp {graph};
    
    size_t choosenVertex {0};
    
   	for (const auto& [i, j] : graph.getAdjacencyList()) {
        sortedVertices.emplace_back(i);
    }

    std::sort(sortedVertices.begin(), sortedVertices.end(),
        [&](size_t a, size_t b) {
            return graph.getVertexDegree(a) > graph.getVertexDegree(b);
    });

    while ((temp.getOrder() > 0) && (choosenVertex < sortedVertices.size())) {
        while (choosenVertex < sortedVertices.size() && 
               !temp.vertexExists(sortedVertices[choosenVertex])) {
            ++choosenVertex;
        }

        if (choosenVertex >= sortedVertices.size()) {
            break;
        }

	    if (temp.getVertexDegree(sortedVertices[choosenVertex]) == 0) {
	        solution.genes[sortedVertices[choosenVertex]] = 3;
	    }
	    
	    else {
	        solution.genes[sortedVertices[choosenVertex]] = 4;

	        for (const auto& it : temp.getAdjacencyList(sortedVertices[choosenVertex])) {
	            if (solution.genes[it] == -1) {
                	solution.genes[it] = 0;
                }
	        }
	    }

	    temp.deleteAdjacencyList(sortedVertices[choosenVertex++]);

	    std::vector<size_t> verticesToRemove;
        
        verticesToRemove.reserve(temp.getOrder());
        
        for (const auto& [i, _]: temp.getAdjacencyList()) {
            if (temp.getVertexDegree(i) == 0) {
                solution.genes[i] = 3;
                verticesToRemove.push_back(i);
            }
        }

        for (const auto& vertex : verticesToRemove) {
            temp.deleteVertex(vertex);
        }
	}
    
    AntColonyOptimization::toggleLabels(graph, solution.genes);
	
    GeneticAlgorithm::fitness(solution);
	
    return solution;
}

bool TripleRomanDomination::feasible(const Graph& graph, const std::vector<int>& solution) {
    bool isValid {false};
    bool hasNeighborAtLeast2 {false};
    
    for (size_t i {0}; i < solution.size(); ++i) {

        isValid = false;

        if (solution[i] == 0) {                                 
            size_t countNeighbors2 {0};
            size_t countNeighbors3 = {0};
            
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
            hasNeighborAtLeast2 = false;
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
