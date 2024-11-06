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

size_t TripleRomanDomination::getGamma3rGeneticAlgorithm() {
    return this->gamma3rGeneticAlgorithm;
}

size_t TripleRomanDomination::getGamma3rACO() {
    return this->gamma3rACO;
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
    Chromosome (*selectedHeuristic)(Graph) = nullptr;

    if (heuristic == 2) 
    	selectedHeuristic = heuristic2;
    else if (heuristic == 3) 
        selectedHeuristic = heuristic3;
    else 
        selectedHeuristic = heuristic1;

    this->geneticAlgorithm->run(geneticAlgorithm->getGenerations(), selectedHeuristic);

    solutionGeneticAlgorithm = this->geneticAlgorithm->getBestSolution();
    std::for_each(solutionGeneticAlgorithm.begin(), solutionGeneticAlgorithm.end(), [&](int element) {
        this->gamma3rGeneticAlgorithm += element;
    });
}

void TripleRomanDomination::runACO() {
   this->ACO->run();
   
   solutionACO = this->ACO->getBestSolution();

   std::for_each(solutionACO.begin(), solutionACO.end(), [&](int element) {
      this->gamma3rACO += element;
   });
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
    std::random_device randomNumber;
    std::mt19937 seed(randomNumber());
    std::uniform_int_distribution<int> gap(0, graph.getOrder() - 1);
	
    size_t choosenVertex = 0;

    while (graph.getOrder() > 0) {
        choosenVertex = gap(seed);
        while (!graph.vertexExists(choosenVertex))
            choosenVertex = gap(seed);
	    
        solution.genes[choosenVertex] = 4;

        for (const auto& it: graph.getAdjacencyList(choosenVertex)) {
            if (solution.genes[it] == -1)
                solution.genes[it] = 0;
        }

        graph.deleteAdjacencyList(choosenVertex);
        if (graph.getOrder() == 1) {
            choosenVertex = graph.getAdjacencyList().begin()->first;
                solution.genes[choosenVertex] = 4;
                graph.deleteVertex(choosenVertex);
        }
    }

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
    std::uniform_int_distribution<int> gap(0, graph.getOrder() - 1);
    
    size_t choosenVertex = 0;
    size_t graphOrder = graph.getOrder();
	
	Graph temp = graph;
	
    while (graph.getOrder() > 0) {
        choosenVertex = gap(seed);
        while (!graph.vertexExists(choosenVertex))
            choosenVertex = gap(seed);

        solution.genes[choosenVertex] = 4;
        for (const auto& it: graph.getAdjacencyList(choosenVertex)) {
            if (solution.genes[it] == -1)
                solution.genes[it] = 0;
        }

        graph.deleteAdjacencyList(choosenVertex);

        for (size_t i = 0; i < graphOrder; ++i) {
            if (graph.vertexExists(i)) {
                if (graph.getVertexDegree(i) == 0) {
                    solution.genes[i] = 3;
                    graph.deleteVertex(i);
                }
            }
        }
    }
    
//    for (size_t i = 0; i < graphOrder; ++i) {
//	    if (solution.genes[i] == 3) {
//	       for (auto& it: temp.getAdjacencyList(i)) 	
//   	          if (solution.genes[it] >= 2) 		              
//          	  		solution.genes[i] = 2; 
//		}
//	}
//
    std::cout << solution << std::endl; 
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

        solution.genes[sortedVertices[choosenVertex]] = 4;

        for (const auto& it : graph.getAdjacencyList(sortedVertices[choosenVertex])) {
            if (solution.genes[it] == -1)
                solution.genes[it] = 0;
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
    
    for (size_t i = 0; i < graphOrder; ++i) {
	    if (solution.genes[i] == 3) {
	       for (auto& it: temp.getAdjacencyList(i)) 	
   	          if (solution.genes[it] >= 2) 		              
          	  		solution.genes[i] = 2; 
		}
	}

    return solution;
}

